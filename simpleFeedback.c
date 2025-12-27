/* simpleFeedback.lv2
 *
 *
 * Copyright (C) 2022 Hans Möller <hmoller@uc.cl>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *
 */

/* include libs */
//#include "stddef.h"
#include "stdint.h"
#include "stdlib.h"
#include "math.h"
#include "lv2.h"
#include "freq_calc.h"

#include "stdio.h"

#define BUFFER_TIME 1 // in sec, multiplied by sample_rates gives buffer_size
#define MIN_FREQ 300
//#define DEBUG  // se supone que lo puse en el makefile y no es necesario

/* class definition */
typedef struct {
    float* in_ptr;
    float* out_ptr;
    float* level_ptr;
    float* attack_ptr;
    float* active_ptr;
    float* harmonic_ptr;

    double rate;  //sample rate

    float* buffer;  // for recording output signal
    float* clean_buffer;  // for recording input signal
    int sample;  // to know where to write in the buffer
    float calc_freq;
    uint8_t prev_active; // to know if it was activated before or not.
    int buffer_size;  // BUFFER_TIME * m->rate
    int delay_pos;

} simpleFeedback;

/* internal core methods */
static LV2_Handle instantiate (const struct LV2_Descriptor *descriptor, double
    sample_rate, const char *bundle_path, const LV2_Feature *const *features){
    simpleFeedback* m = (simpleFeedback*)calloc(1, sizeof(simpleFeedback));
    if (!m) {
        return NULL;
    }

    m->rate        = sample_rate;
    m->buffer_size = (uint32_t)(BUFFER_TIME * sample_rate);

    if (m->buffer_size == 0) {
        free(m);
        return NULL;
    }

    m->buffer = (float*)calloc(m->buffer_size, sizeof(float));
    m->clean_buffer = (float*)calloc(m->buffer_size, sizeof(float));

    if (!m->buffer || !m->clean_buffer) {
        free(m->buffer);
        free(m->clean_buffer);
        free(m);
        return NULL;
    }

    m->sample     = 0;
    m->calc_freq  = 0.0f;
    m->delay_pos  = 0;
    m->prev_active = 0;

    return m;
}


static void connect_port (LV2_Handle instance, uint32_t port, void 
    *data_location){
    simpleFeedback* m = (simpleFeedback*) instance;
    if (!m) return;

    switch (port){
    case 0:
        m->in_ptr = (float*) data_location;
        break;
    case 1:
        m->out_ptr = (float*) data_location;
        break;
    case 2:
        m->level_ptr = (float*) data_location;
        break;
    case 3:
        m->active_ptr = (float*) data_location;
        break;
    case 4:
        m->attack_ptr = (float*) data_location;
        break;
    case 5:
        m->harmonic_ptr = (float*) data_location;
        break;
    default:
        break;
    }
}

static void activate (LV2_Handle instance){
simpleFeedback* m = (simpleFeedback*)instance;
    if (!m) return;
    
    /*for (int i = 0; i < m->buffer_size; i++) {
		m->buffer[i] = 0;
        m->clean_buffer[i] = 0;
	}*/
    if (m->buffer && m->clean_buffer && m->buffer_size > 0) {
        memset(m->buffer,       0, m->buffer_size * sizeof(float));
        memset(m->clean_buffer, 0, m->buffer_size * sizeof(float));
    }

    m->sample    = 0;
    m->calc_freq = 0.0f;
    m->delay_pos = 0;

    if (m->active_ptr && *m->active_ptr < 0.5f)
        m->prev_active = 0;
    else
        m->prev_active = 1;

    /* Para evitar warning si attack_ptr aún no se usa */
    (void)m->attack_ptr;
}

static void run(LV2_Handle instance, uint32_t sample_count)
{
    simpleFeedback* m = (simpleFeedback*)instance;
    if (!m) return;

    if ((!m->in_ptr) || (!m->out_ptr) || (!m->level_ptr) ||
        (!m->active_ptr) || (!m->harmonic_ptr) || (!m->attack_ptr) ||
        (!m->buffer) || (!m->clean_buffer) || (m->buffer_size == 0)) {
        return;
    }

    float temp_freq = 0.0f;
    const uint32_t buf_size = m->buffer_size;  // dado que se ocupa muchas veces es más rápido leer la memoria que el puntero. Además al saber el compilador q la variable no cambia puede introducir optimizaciones.

    if (*m->active_ptr < 0.5f) {
        /* Bypass + guardado de clean_buffer */
        for (uint32_t i = 0; i < sample_count; ++i) {
            m->out_ptr[i] = m->in_ptr[i];
            m->clean_buffer[m->sample] = m->in_ptr[i];
            m->sample = (m->sample + 1U) % buf_size;  // como m->sample es uint, le ponemos U al 1 para que sea unsigned.
        }

        if (m->prev_active == 1) {
            /* Si antes estaba activo, vaciar buffer y resetear índice */
            for (uint32_t i = 0; i < buf_size; ++i) {
                m->buffer[i] = 0.0f;
            }
            m->sample = 0;
            m->prev_active = 0;
        }
        return;
    }

    /* Estado activo */
    if (m->prev_active == 0) {
        m->prev_active = 1;
    }

    /* Solo recalculamos frecuencia si tenemos una ventana mínima */
    const uint32_t min_samples = (uint32_t)(4.0 * m->rate / (double)MIN_FREQ);
    if (m->sample > min_samples) {
        temp_freq = (float)fft_autocorr_freq(
            m->clean_buffer, m->sample, buf_size, m->rate);

        if (temp_freq > 20.0f) {
            m->calc_freq = temp_freq;

            /* Proteger harmonic_ptr por si el host manda 0 */
            float harmonic = (*m->harmonic_ptr > 0.01f) ? *m->harmonic_ptr : 1.0f;

            double delay = m->rate / (double)m->calc_freq / (double)harmonic;
            if (delay < 1.0) {
                delay = 1.0;
            }
            if (delay > buf_size - 1) {
                delay = buf_size - 1;
            }

            m->delay_pos = (uint32_t)delay;

            /* No hacer printf en tiempo real */
            /* printf("%f\n", m->calc_freq); */
            #ifdef DEBUG
            log_msg(m, m->log_Notice, "freq=%f\n", m->calc_freq);
            #endif
        }
    }

    for (uint32_t i = 0; i < sample_count; ++i) {
        //calculate which position we must read from buffer
        uint32_t eco_pos = (m->sample + buf_size - m->delay_pos) % buf_size;

        float delayed = m->buffer[eco_pos];
        float in      = m->in_ptr[i];

        m->out_ptr[i] = in + delayed * (*m->level_ptr);

        /* Guardamos señal con feedback en buffer, y la señal limpia en clean_buffer */
        m->buffer[m->sample]       = in + delayed;
        m->clean_buffer[m->sample] = in;

        m->sample = (m->sample + 1U) % buf_size;
    }
}

static void deactivate(LV2_Handle instance)
{
    /* En este diseño no liberamos nada aquí
       (los buffers se mantienen hasta cleanup).
       Importante: no hacer free en un posible contexto RT. */
    (void)instance;
}

static void cleanup(LV2_Handle instance)
{
    simpleFeedback* m = (simpleFeedback*)instance;
    if (!m) return;

    if (m->buffer) {
        free(m->buffer);
        m->buffer = NULL;
    }
    if (m->clean_buffer) {
        free(m->clean_buffer);
        m->clean_buffer = NULL;
    }

    free(m);
}

static const void* extension_data(const char* uri)
{
    (void)uri;
    return NULL;
}

/* descriptor */
static const LV2_Descriptor descriptor = {
    "https://github.com/hmollercl/simpleFeedback",
    instantiate,
    connect_port,
    activate,
    run,
    deactivate,
    cleanup,
    extension_data
};

LV2_SYMBOL_EXPORT const LV2_Descriptor* lv2_descriptor(uint32_t index)
{
    if (index == 0) {
        return &descriptor;
    }
    return NULL;
}