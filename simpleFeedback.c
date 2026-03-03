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
#include "string.h"
#include "math.h"
#include "lv2.h"
//#include "freq_calc.h"

#include "math.h"
#include <sys/types.h>
#include <float.h>
#include <fftw3.h>
#define F_MIN 100.0f     // Hz (Mi2 ~82 Hz, nota más baja)
#define F_MAX 1400.0f    // Hz (trastes altos, agudos)
#define FFT_SIZE 16384   // tamaño de la ventana para calcular fft


#define BUFFER_TIME 1 // in sec, multiplied by sample_rates gives buffer_size
#define MIN_FREQ 300
//#define DEBUG  // in makefile

/* class definition */
typedef struct {
    float* in_ptr;
    float* out_ptr;
    float* level_ptr;
    float* attack_ptr;
    float* active_ptr;
    float* harmonic_ptr;
    float* lock_ptr;

    double rate;  //sample rate

    float* buffer;  // for recording output signal
    float* clean_buffer;  // for recording input signal
    int sample;  // to know where to write in the buffer
    float calc_freq;
    uint8_t prev_active; // to know if it was activated before or not.
    int buffer_size;  // BUFFER_TIME * m->rate
    float delay_pos;
    float delay_prev_pos;
    float delay_xfade;
    float delay_xfade_inc;
    uint8_t delay_xfade_active;
    float effect_gain;   // ganancia “suavizada” del efecto, para usar attack_ptr
    uint8_t freq_locked;
    uint8_t onset_flag;
    float onset_env;
    float bp_b0;
    float bp_b1;
    float bp_b2;
    float bp_a1;
    float bp_a2;
    float bp_z1;
    float bp_z2;
    float bp_center_hz;

    // === Pitch detection (autocorrelación) ===
    float*          pd_time;      // buffer de entrada para la FFT (FFT_SIZE)
    fftwf_complex*  pd_freq;      // salida FFT compleja (FFT_SIZE/2+1)
    float*          pd_autocorr;  // autocorrelación (FFT_SIZE)
    fftwf_plan      pd_plan_fwd;  // plan FFT r2c
    fftwf_plan      pd_plan_inv;  // plan IFFT c2r

} simpleFeedback;

static inline float
clamp_delay_pos(float delay_pos, uint32_t buf_size)
{
    if (delay_pos < 1.0f) {
        return 1.0f;
    }
    if (delay_pos > (float)(buf_size - 1U)) {
        return (float)(buf_size - 1U);
    }
    return delay_pos;
}

static inline float
read_delay_linear(const float* buffer, uint32_t write_pos, uint32_t buf_size, float delay_pos)
{
    float read_pos = (float)write_pos - delay_pos;
    if (read_pos < 0.0f) {
        read_pos += (float)buf_size;
    }

    uint32_t pos0 = (uint32_t)read_pos;
    uint32_t pos1 = (pos0 + 1U) % buf_size;
    float frac = read_pos - (float)pos0;
    return buffer[pos0] + frac * (buffer[pos1] - buffer[pos0]);
}

static inline float
soft_clip_fast(float x)
{
    float y = 1.4f * x;
    return y / (1.0f + fabsf(y));
}

static inline void
update_bandpass_coeffs(simpleFeedback* m, float center_hz)
{
    const float q = 8.0f;
    const float min_hz = 20.0f;
    const float max_hz = 0.45f * (float)m->rate;

    if (center_hz < min_hz) {
        center_hz = min_hz;
    } else if (center_hz > max_hz) {
        center_hz = max_hz;
    }

    float w0 = 2.0f * (float)M_PI * center_hz / (float)m->rate;
    float cos_w0 = cosf(w0);
    float sin_w0 = sinf(w0);
    float alpha = sin_w0 / (2.0f * q);

    float a0 = 1.0f + alpha;
    m->bp_b0 = alpha / a0;
    m->bp_b1 = 0.0f;
    m->bp_b2 = -alpha / a0;
    m->bp_a1 = (-2.0f * cos_w0) / a0;
    m->bp_a2 = (1.0f - alpha) / a0;
    m->bp_center_hz = center_hz;
}

static inline float
process_bandpass(simpleFeedback* m, float x)
{
    float y = m->bp_b0 * x + m->bp_z1;
    m->bp_z1 = m->bp_b1 * x - m->bp_a1 * y + m->bp_z2;
    m->bp_z2 = m->bp_b2 * x - m->bp_a2 * y;
    return y;
}

static float
autocorr_freq_rt(simpleFeedback* m,
                 const float*    buffer,
                 int             position,
                 int             buffer_size)
{
    if (!m || !buffer || buffer_size <= 0) {
        return 0.0f;
    }

    const int N = FFT_SIZE;

    // === 1. Copiar N muestras del buffer circular a m->pd_time ===
    // Suponemos que 'position' es la próxima posición de escritura,
    // y que los datos válidos son los N samples anteriores.
    int start = position - N;
    while (start < 0) {
        start += buffer_size;
    }

    int idx = start;
    float sumsq = 0.0f;  //for RMS calc
    for (int i = 0; i < N; ++i) {
        m->pd_time[i] = buffer[idx];
        float x = buffer[idx];
        sumsq += x * x;
        idx++;
        if (idx >= buffer_size) {
            idx = 0;
        }
    }

    float rms = sqrtf(sumsq / (float)N);

    // mín RMS (adjustable)
    const float RMS_MIN = 0.001f;   // to use between 0.001–0.005
    if (rms < RMS_MIN) {
        return 0.0f;
    }

    // === 2. FFT (forward) ===
    fftwf_execute(m->pd_plan_fwd);

    // === 3. Espectro de potencia: |X(f)|^2 ===
    int spec_size = N / 2 + 1;
    for (int i = 0; i < spec_size; ++i) {
        float re = m->pd_freq[i][0];
        float im = m->pd_freq[i][1];
        m->pd_freq[i][0] = re * re + im * im; // real power
        m->pd_freq[i][1] = 0.0f;              // imag power = 0
    }

    // === 4. IFFT → autocorrelation ===
    fftwf_execute(m->pd_plan_inv);

    // === 5. Normalizar autocorrelación ===
    for (int i = 0; i < N; ++i) {
        m->pd_autocorr[i] /= (float)N;
    }

    // Energía en tau = 0
    float r0 = m->pd_autocorr[0];
    if (r0 <= 1e-9f) {
        return 0.0f;
    }

    // === 6. Buscar lag en rango [lag_min, lag_max] ===
    int lag_min = (int)(m->rate / F_MAX);  // lag mínimo (freq más alta)
    int lag_max = (int)(m->rate / F_MIN);  // lag máximo (freq más baja)

    if (lag_min < 1) lag_min = 1;
    if (lag_max >= N) lag_max = N - 1;
    if (lag_min >= lag_max) {
        return 0.0f;
    }

    float max_val = -1.0e30f;
    int   max_lag = lag_min;

    for (int lag = lag_min; lag <= lag_max; ++lag) {
        float v = m->pd_autocorr[lag];
        if (v > max_val) {
            max_val = v;
            max_lag = lag;
        }
    }

    if (max_lag <= 0) {
        return 0.0f;
    }

    // Si el mejor pico cae en el borde del rango buscado, suele ser una
    // detección ambigua/saturada en F_MIN o F_MAX. En ese caso descartamos
    // la estimación para mantener la frecuencia previa en run().
    if (max_lag == lag_min || max_lag == lag_max) {
        return 0.0f;
    }

    // === 7. Evaluar calidad del pico (normalizado) ===
    float norm_peak = max_val / r0;

    // Umbral de "periodicidad" (ajustable: 0.2–0.4 suelen ir bien)
    const float CORR_MIN = 0.25f;

    if (norm_peak < CORR_MIN) {
        // Pico de autocorrelación débil → probablemente ruido / sin pitch claro
        return 0.0f;
    }

    // === 8. Lag → frecuencia ===
    float freq = (float)m->rate / (float)max_lag;
    return freq;


    /*Cómo ajustar los umbrales en la práctica
    RMS_MIN:
    Empieza con algo bajo, tipo 0.001 o 0.002.
    Si sigue actualizando pitch en colas muy débiles → súbelo un poco.

    CORR_MIN:
    Empieza con 0.25.
    Si te da muchos falsos pitches con ruido → sube a 0.3–0.35.
    Si se vuelve demasiado “perezoso” para detectar notas suaves → bájalo a 0.2.*/
}


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

    // ==== Pitch detection: mem for FFTW ====
    m->pd_time = (float*)fftwf_malloc(sizeof(float) * FFT_SIZE);
    m->pd_freq = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (FFT_SIZE/2 + 1));
    m->pd_autocorr = (float*)fftwf_malloc(sizeof(float) * FFT_SIZE);

    if (!m->pd_time || !m->pd_freq || !m->pd_autocorr) {
        if (m->pd_time) fftwf_free(m->pd_time);
        if (m->pd_freq) fftwf_free(m->pd_freq);
        if (m->pd_autocorr) fftwf_free(m->pd_autocorr);
        free(m->buffer);
        free(m->clean_buffer);
        free(m);
        return NULL;
    }

    // ==== Create FFTW plan (no-RT, here MEASURE can be used) ====
    m->pd_plan_fwd = fftwf_plan_dft_r2c_1d(
        FFT_SIZE, m->pd_time, m->pd_freq, FFTW_MEASURE);
    m->pd_plan_inv = fftwf_plan_dft_c2r_1d(
        FFT_SIZE, m->pd_freq, m->pd_autocorr, FFTW_MEASURE);

    if (!m->pd_plan_fwd || !m->pd_plan_inv) {
        if (m->pd_plan_fwd) fftwf_destroy_plan(m->pd_plan_fwd);
        if (m->pd_plan_inv) fftwf_destroy_plan(m->pd_plan_inv);
        fftwf_free(m->pd_time);
        fftwf_free(m->pd_freq);
        fftwf_free(m->pd_autocorr);
        free(m->buffer);
        free(m->clean_buffer);
        free(m);
        return NULL;
    }

    m->sample     = 0;
    m->calc_freq  = 0.0f;
    m->delay_pos  = 1.0f;
    m->delay_prev_pos = 1.0f;
    m->delay_xfade = 1.0f;
    m->delay_xfade_inc = 0.0f;
    m->delay_xfade_active = 0;
    m->freq_locked = 0;
    m->onset_flag = 0;
    m->onset_env = 0.0f;
    m->bp_z1      = 0.0f;
    m->bp_z2      = 0.0f;
    update_bandpass_coeffs(m, 1000.0f);
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
    case 6:
        m->lock_ptr = (float*) data_location;
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
    m->delay_pos = 1.0f;
    m->delay_prev_pos = 1.0f;
    m->delay_xfade = 1.0f;
    m->delay_xfade_inc = 0.0f;
    m->delay_xfade_active = 0;
    m->effect_gain = 0.0f;
    m->freq_locked = 0;
    m->onset_flag = 0;
    m->onset_env = 0.0f;
    m->bp_z1 = 0.0f;
    m->bp_z2 = 0.0f;

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
            m->delay_xfade_active = 0;
            m->delay_xfade = 1.0f;
            m->delay_xfade_inc = 0.0f;
            m->freq_locked = 0;
            m->onset_flag = 0;
            m->onset_env = 0.0f;
            m->bp_z1 = 0.0f;
            m->bp_z2 = 0.0f;
            m->prev_active = 0;
        }
        return;
    }

    /* Estado activo */
    if (m->prev_active == 0) {
        m->prev_active = 1;
        m->effect_gain = 0.0f;  // when activate start with 0.
        m->delay_xfade_active = 0;
        m->delay_xfade = 1.0f;
    }

    /* Solo recalculamos frecuencia si tenemos una ventana mínima */
    const uint32_t min_samples = (uint32_t)(4.0 * m->rate / (double)MIN_FREQ);
    const float xfade_time_sec = 0.02f;
    const float lock_in_rel = 0.06f;
    const float lock_out_rel = 0.12f;
    uint8_t onset_now = m->onset_flag;
    m->onset_flag = 0;

    if (m->sample > min_samples) {
        temp_freq = (float)autocorr_freq_rt(m,
                                       m->clean_buffer,
                                       m->sample,
                                       m->buffer_size);

        if (temp_freq > 0.0f) {
            uint8_t accept_update = 0;
            uint8_t lock_enabled = (!m->lock_ptr || *m->lock_ptr >= 0.5f);
            if (!lock_enabled) {
                accept_update = 1;
                m->freq_locked = 0;
            } else if (m->calc_freq <= 0.0f) {
                accept_update = 1;
                m->freq_locked = 1;
            } else {
                float rel_diff = fabsf(temp_freq - m->calc_freq) / fmaxf(m->calc_freq, 1.0f);
                if (m->freq_locked) {
                    if (rel_diff <= lock_out_rel || onset_now) {
                        accept_update = 1;
                    } else {
                        m->freq_locked = 0;
                    }
                } else if (rel_diff <= lock_in_rel || onset_now) {
                    accept_update = 1;
                    m->freq_locked = 1;
                }
            }

            if (accept_update) {
                // Hay un pitch confiable y estable → actualiza con suavizado
                const float alpha = 0.3f; // 0 = ultra suave, 1 = sin suavizado
                if (m->calc_freq <= 0.0f) {
                    m->calc_freq = temp_freq;
                } else {
                    m->calc_freq = alpha * temp_freq + (1.0f - alpha) * m->calc_freq;
                }

                float harmonic = (*m->harmonic_ptr > 0.01f ? *m->harmonic_ptr : 1.0f);
                float target_hz = m->calc_freq * harmonic;
                if (target_hz > 1.0f) {
                    float new_delay_pos = clamp_delay_pos((float)(m->rate / target_hz), buf_size);
                    if (fabsf(new_delay_pos - m->delay_pos) > 0.0001f) {
                        m->delay_prev_pos = m->delay_pos;
                        m->delay_pos = new_delay_pos;
                        m->delay_xfade = 0.0f;
                        m->delay_xfade_inc = 1.0f / fmaxf(1.0f, xfade_time_sec * (float)m->rate);
                        m->delay_xfade_active = 1;
                    }
                }

                if (fabsf(target_hz - m->bp_center_hz) > 0.5f) {
                    update_bandpass_coeffs(m, target_hz);
                }

                /* No hacer printf en tiempo real */
                /* printf("%f\n", m->calc_freq); */
                #ifdef DEBUG
                printf("%f\n", m->calc_freq);
                #endif
            }
        }
    }

    // ==== ATTACK (fade-in del nivel de efecto) ====
    float target_level = *m->level_ptr;
    if (target_level < 0.0f) target_level = 0.0f;

    float attack_norm = *m->attack_ptr;
    if (attack_norm < 0.0f) attack_norm = 0.0f;
    if (attack_norm > 1.0f) attack_norm = 1.0f;

    const float ATK_MIN = 0.02f;
    const float ATK_MAX = 2.0f;
    float attack_time = ATK_MIN * powf(ATK_MAX / ATK_MIN, attack_norm);

    float step;
    if (attack_time <= 0.000001f) {
        step = target_level;
    } else {
        step = target_level / (attack_time * (float)m->rate);
        if (step > target_level) step = target_level;
    }

    float block_peak = 0.0f;
    for (uint32_t i = 0; i < sample_count; ++i) {

        // ramp lineal hacia target_level
        if (m->effect_gain < target_level) {
            m->effect_gain += step;
            if (m->effect_gain > target_level)
                m->effect_gain = target_level;
        } else {
            m->effect_gain = target_level;
        }

        float delayed = read_delay_linear(m->buffer, (uint32_t)m->sample, buf_size, m->delay_pos);
        if (m->delay_xfade_active) {
            float delayed_prev = read_delay_linear(m->buffer, (uint32_t)m->sample, buf_size, m->delay_prev_pos);
            delayed = delayed_prev + (delayed - delayed_prev) * m->delay_xfade;
            m->delay_xfade += m->delay_xfade_inc;
            if (m->delay_xfade >= 1.0f) {
                m->delay_xfade = 1.0f;
                m->delay_xfade_active = 0;
                m->delay_prev_pos = m->delay_pos;
            }
        }

        float in      = m->in_ptr[i];
        float abs_in = fabsf(in);
        if (abs_in > block_peak) {
            block_peak = abs_in;
        }

        m->out_ptr[i] = in + delayed * m->effect_gain;

        // feedback interno también con attack (clave)
        float feedback_in = in + delayed * m->effect_gain;
        float clipped = soft_clip_fast(feedback_in);
        m->buffer[m->sample] = process_bandpass(m, clipped);

        m->clean_buffer[m->sample] = in;
        m->sample = (m->sample + 1U) % buf_size;
    }

    // onset simple por bloque (para permitir relock cuando hay ataque nuevo)
    if (block_peak > 0.03f && block_peak > 1.8f * m->onset_env) {
        m->onset_flag = 1;
    }
    m->onset_env = 0.995f * m->onset_env + 0.005f * block_peak;

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

    // Pitch detection (FFTW)
    if (m->pd_plan_fwd) {
        fftwf_destroy_plan(m->pd_plan_fwd);
        m->pd_plan_fwd = NULL;
    }
    if (m->pd_plan_inv) {
        fftwf_destroy_plan(m->pd_plan_inv);
        m->pd_plan_inv = NULL;
    }

    if (m->pd_time) {
        fftwf_free(m->pd_time);
        m->pd_time = NULL;
    }
    if (m->pd_freq) {
        fftwf_free(m->pd_freq);
        m->pd_freq = NULL;
    }
    if (m->pd_autocorr) {
        fftwf_free(m->pd_autocorr);
        m->pd_autocorr = NULL;
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
