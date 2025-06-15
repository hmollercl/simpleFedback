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

#define BUFFER_SIZE 3 // multiplied by sample_rates, 3 -> 3sec.

/* class definition */
typedef struct {
    float* in_ptr;
    float* out_ptr;
    float* level_ptr;
    float* delay_ptr;
    float* attack_ptr;
    float* active_ptr;
    float* freq_ptr;
    float* q_ptr;

    double rate;  //sample rate

    float* buffer;  // for recording
    int sample;  // to know where to read the buffer

    float x[3];  //filter values
	float y[3];  //filter values

} simpleFeedback;

/* internal core methods */
static LV2_Handle instantiate (const struct LV2_Descriptor *descriptor, double
    sample_rate, const char *bundle_path, const LV2_Feature *const *features){
    simpleFeedback* m = (simpleFeedback*) calloc (1, sizeof (simpleFeedback));
    if(m)
        m->rate = sample_rate;
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
        m->delay_ptr = (float*) data_location;
        break;
    case 4:
        m->active_ptr = (float*) data_location;
        break;
    case 5:
        m->attack_ptr = (float*) data_location;
        break;
    case 6:
        m->freq_ptr = (float*) data_location;
        break;
    case 7:
        m->q_ptr = (float*) data_location;
        break;
    default:
        break;
    }
}

static void activate (LV2_Handle instance){
    simpleFeedback* m = (simpleFeedback*) instance;

	m->buffer = malloc(m->rate * BUFFER_SIZE * sizeof(float));
	for (int i = 0; i < m->rate * BUFFER_SIZE; i++) {
		m->buffer[i] = 0;
	}
	m->sample = 0;

    for (int i = 0; i < 3; i++) {
		m->x[i] = 0;
		m->y[i] = 0;
	}
}

static void run (LV2_Handle instance, uint32_t sample_count){
    simpleFeedback* m = (simpleFeedback*) instance;
    //const float        freq      = 300;
	//const float        q         = 1;

    //banpdass biquad filter
	float w0 = 2 * 3.1416 * *m->freq_ptr / m->rate;  // TODO calculate freq
	float alpha = sin(w0) / (2 * *m->q_ptr);  // TODO define q
	float b0 = (1 - cos(w0)) / 2;
	float b1 = 1 - cos(w0);
	float b2 = (1 - cos(w0)) / 2;
	float a0 = 1 + alpha;
	float a1 = -2 * cos(w0);
	float a2 = 1 - alpha;

    if (!m) return;
    if ((!m->in_ptr) || (!m->out_ptr) || (!m->level_ptr) || 
        (!m->active_ptr) || (!m->delay_ptr) || (!m->attack_ptr)) return;

    if (*m->active_ptr < 0.5) { // or active_state == false if using boolean
        // Bypass: Copy input to output
        for (uint32_t i = 0; i < sample_count; ++i) {
            m->out_ptr[i] = m->in_ptr[i];
        }
    }
    else{
        uint32_t eco_pos;
        for (uint32_t i = 0; i < sample_count; i++) {
            //calculate which position we must read from buffer.
            if (m->sample < (*m->delay_ptr * m->rate))
                eco_pos = m->sample + m->rate * BUFFER_SIZE - *m->delay_ptr * m->rate;
            else
                eco_pos = (uint32_t)(m->sample - (*m->delay_ptr * m->rate));

            // output = input + eco
            //m->out_ptr[i] = m->in_ptr[i] + m->buffer[eco_pos] * *m->level_ptr;
            
            //aplico el filtro
            m->x[2] = m->x[1]; // x [z-2]
            m->x[1] = m->x[0]; // x [z-1]
            //m->x[0] = input[pos]; // x [z]
            m->x[0] = m->in_ptr[i] + m->buffer[eco_pos] * *m->level_ptr;
            m->y[2] = m->y[1]; // y [z-2]
            m->y[1] = m->y[0]; // y [z-1]
            m->y[0] = (b0/a0) * m->x[0] + (b1/a0) * m->x[1] + (b2/a0) * m->x[2]
                                                    - (a1/a0) * m->y[1] - (a2/a0) * m->y[2];
            //output[pos] = m->y[0];
            m->out_ptr[i] = m->y[0];

            // save output in buffer
            m->buffer[m->sample] = m->out_ptr[i];

            //calculate next sample where we will write in buffer.
            m->sample++;
            if (m->sample > m->rate * BUFFER_SIZE)
                m->sample = 0;
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
