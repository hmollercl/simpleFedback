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

    float x[3];  //filter values  // not necessary
	float y[3];  //filter values  // not necessary

} simpleFeedback;

/* internal core methods */
static LV2_Handle instantiate (const struct LV2_Descriptor *descriptor, double
    sample_rate, const char *bundle_path, const LV2_Feature *const *features){
    simpleFeedback* m = (simpleFeedback*) calloc (1, sizeof (simpleFeedback));
    if(m)
        m->rate = sample_rate;
        m->buffer_size = BUFFER_TIME * m->rate;
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
    simpleFeedback* m = (simpleFeedback*) instance;

	m->buffer = malloc(m->buffer_size * sizeof(float));
    m->clean_buffer = malloc(m->buffer_size * sizeof(float));
	for (int i = 0; i < m->buffer_size; i++) {
		m->buffer[i] = 0;
        m->clean_buffer[i] = 0;
	}
	m->sample = 0;
    m->calc_freq = 0;
    m->delay_pos=0;

    if(*m->active_ptr < 0.5)
        m->prev_active = 0;
    else
        m->prev_active = 1;

    for (int i = 0; i < 3; i++) {
		m->x[i] = 0;
		m->y[i] = 0;
	}
}

static void run (LV2_Handle instance, uint32_t sample_count){
    simpleFeedback* m = (simpleFeedback*) instance;
    float temp_freq = 0;
    
    if (!m) return;
    if ((!m->in_ptr) || (!m->out_ptr) || (!m->level_ptr) || 
        (!m->active_ptr) || (!m->harmonic_ptr) || (!m->attack_ptr)) return;

    //TODO creo que acá es más simple si no está active escribir 0s.
    if (*m->active_ptr < 0.5) { // or active_state == false if using boolean
        // Bypass: Copy input to output
        for (uint32_t i = 0; i < sample_count; ++i) {
            m->out_ptr[i] = m->in_ptr[i];
            m->clean_buffer[m->sample] = m->in_ptr[i];
            m->sample = (m->sample + 1) % (int)(m->buffer_size);
        }
        if (m->prev_active == 1){
            //if before was activated, clear buffer and restart sample position
            for (int i = 0; i < m->buffer_size; i++) {
                m->buffer[i] = 0;
            }
            m->sample = 0;
            m->prev_active = 0;
        }
    }
    else{
        if (m->prev_active == 0)
            m->prev_active = 1;
        uint32_t eco_pos;
        if (m->sample > (4 * m->rate / 300)){ //because min_win_length is rate/300 *2 for 300Hz min freq.
            //temp_freq = (float) (fft_freq(m->clean_buffer, m->sample, m->buffer_size, m->rate));
            temp_freq = (float) (fft_autocorr_freq(m->clean_buffer, m->sample, m->buffer_size, m->rate));
            //temp_freq = (float) (yin_freq(m->clean_buffer, m->sample, m->buffer_size, m->rate));
            printf("%f \n",temp_freq);
            
            if (temp_freq > 20){
                m->calc_freq = temp_freq;
                m->delay_pos = m->rate / m->calc_freq / *m->harmonic_ptr;  // to match delay with frequency

                //printf("%f \n",m->calc_freq);
            }
        }

        for (uint32_t i = 0; i < sample_count; i++) {
            //calculate which position we must read from buffer
            eco_pos = (m->sample - m->delay_pos + m->buffer_size) % m->buffer_size;
    
            m->out_ptr[i] = m->in_ptr[i] + m->buffer[eco_pos] * *m->level_ptr;

            // save output in buffer
            m->buffer[m->sample] = m->out_ptr[i];
            m->clean_buffer[m->sample] = m->in_ptr[i];

            //calculate next sample where we will write in buffer.
            /*m->sample++;
            if (m->sample > m->rate * BUFFER_SIZE)
                m->sample = 0;*/
            m->sample = (m->sample + 1) % (int)(m->buffer_size);
        }
    }
}

static void deactivate (LV2_Handle instance)
{
    /* not needed here */
}

static void cleanup (LV2_Handle instance)
{
    simpleFeedback* m = (simpleFeedback*) instance;
    if (!m) return;
    free (m);
}

static const void* extension_data (const char *uri)
{
    return NULL;
}

/* descriptor */
static LV2_Descriptor const descriptor =
{
    "https://github.com/hmollercl/simpleFeedback",
    instantiate,
    connect_port,
    activate /* or NULL */,
    run,
    deactivate /* or NULL */,
    cleanup,
    extension_data /* or NULL */
};

/* interface */
const LV2_SYMBOL_EXPORT LV2_Descriptor* lv2_descriptor (uint32_t index)
{
    if (index == 0) return &descriptor;
    else return NULL;
}
