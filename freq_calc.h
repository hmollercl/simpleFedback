#include "math.h"
#include <sys/types.h>
#include <float.h>

#include <fftw3.h>

float fft_freq(float *buffer, int position, int buffer_size, int rate){
    // quizá se debiera reordenar el buffer...
    // Crear el plan FFTW
    fftwf_complex *out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * (buffer_size/2 + 1));
    fftwf_plan p = fftwf_plan_dft_r2c_1d(buffer_size, buffer, out, FFTW_ESTIMATE);
    
    // Realizar la FFT
    fftwf_execute(p);   
    // Encontrar la frecuencia máxima
    double max_amplitude = 0;
    int max_index = 0;
    
    for (int i = 0; i < buffer_size/2 + 1; i++) {
        double amplitude = sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]);
        if (amplitude > max_amplitude) {
            max_amplitude = amplitude;
            max_index = i;
        }
    }
    
    // Liberar recursos
    fftwf_destroy_plan(p);
    fftwf_free(out);
    return(max_index * rate / buffer_size);
}



float asdf_freq(float* buffer, int position, int buffer_size, int rate){
    //AMDF algorithm for pitch detection
    float sensitivity = 1;
    int window_length = (int) 1.9 * rate / 300; // 300Hz min freq to look for.
    //TODO organize buffer so position is the last
    //for now, we are asumming that position is greater than 2 x window_length
    // asuming 48k sampling rate in some calculations

    float best_result = 9999999;
    u_int8_t best_diff;
    float temp_result;
    u_int8_t i;
    u_int8_t j;

    //TODO rewrite this with only buffer and window_length and calculate start i
    // i = (int) rate / 1600; 
    for (i = 30; i <= window_length; i++){ //start with 30 for max ~1.6kHz, (sample_rate/i max freq),
        temp_result = 0;
        for (j = 0; j <= window_length; j++){
            //temp_result += fabs(buffer[position-1-j] - buffer[position-1-i-j]);
            temp_result += pow(buffer[position-j] - buffer[position-i-j], 2);
        }
        temp_result = temp_result/(window_length * pow(0.99, window_length - i));
        if (temp_result < best_result)
        {
            best_result = temp_result;
            best_diff = i;
        }
    }
    if (best_result < sensitivity)
        return rate/best_diff;
    else
        return 99999;
}

float amdf_freq(float* buffer, int position, int buffer_size, int rate){
    //AMDF algorithm for pitch detection
    float sensitivity = 1;
    int window_length = (int) 1.9 * rate / 300; // 300Hz min freq to look for.
    //TODO organize buffer so position is the last
    //for now, we are asumming that position is greater than 2 x window_length
    // asuming 48k sampling rate in some calculations

    float best_result = 9999999;
    u_int8_t best_diff;
    float temp_result;
    u_int8_t i;
    u_int8_t j;

    //TODO rewrite this with only buffer and window_length and calculate start i
    // i = (int) rate / 1600; 
    for (i = 30; i <= window_length; i++){ //start with 30 for max ~1.6kHz, (sample_rate/i max freq),
        temp_result = 0;
        for (j = 0; j < window_length; j++){
            //temp_result += fabs(buffer[position-1-j] - buffer[position-1-i-j]);
            temp_result += fabs(buffer[position-j] - buffer[position-i-j]);
        }
        //temp_result = temp_result/window_length;
        temp_result = temp_result/(window_length * pow(0.99, -window_length + i));
        if (temp_result < best_result)
        {
            best_result = temp_result;
            best_diff = i;
        }
    }
    if (best_result < sensitivity)
        return rate/best_diff;
    else
        return 0;
}

float acf_freq(float* buffer, int position, int buffer_size, int rate){
    //ACF algorithm for pitch detection
    float sensitivity = 0.01;
    int window_length = (int) 1.9 * rate / 300; // 300Hz min freq to look for.
    //TODO organize buffer so position is the last
    //for now, we are asumming that position is greater than 2 x window_length
    // asuming 48k sampling rate in some calculations

    float best_result = 0;
    u_int8_t best_diff;
    float temp_result;
    u_int8_t i;
    u_int8_t j;

    //TODO rewrite this with only buffer and window_length and calculate start i
    // i = (int) rate / 1600; 
    for (i = 30; i <= window_length; i++){ //start with 30 for max ~1.6kHz, (sample_rate/i max freq),
        temp_result = 0;
        for (j = 0; j <= window_length; j++){
            temp_result += buffer[position-j] * buffer[position-i-j];
        }
        temp_result = temp_result/window_length;
        if (temp_result >= best_result)
        {
            best_result = temp_result;
            best_diff = i;
        }
    }
    return rate/best_diff;
}



#ifndef YIN_THRESHOLD
#define YIN_THRESHOLD 0.10f     // Umbral típico: 0.10–0.15
#endif

#ifndef YIN_MIN_FREQUENCY
#define YIN_MIN_FREQUENCY 300.0f // Para guitarra puedes usar 82.0f (E2)
#endif

#ifndef YIN_MAX_FREQUENCY
#define YIN_MAX_FREQUENCY 2000.0f
#endif
// -------------------------------------------------------------

static inline float yin_parabolic_min_pos(float ym1, float y0, float yp1) {
    // Devuelve offset submuestral del mínimo (≈ -0.5..0.5)
    float denom = (ym1 - 2.0f*y0 + yp1);
    if (fabsf(denom) < 1e-12f) return 0.0f;
    return 0.5f * (ym1 - yp1) / denom;
}

// Firma solicitada:
float yin_freq(float* buffer, int position, int buffer_size, int rate) {
    if (!buffer || buffer_size < 4 || rate <= 0) return 0.0f;

    //const float *x = buffer + position;      // comienzo de la ventana
    //revisar esto...
    //const int N = buffer_size;               // largo de la ventana
    //const float sr = (float)rate;
    const int N = (int) 1.9 * rate / 300; // 300Hz min freq to look for.

    // Rango de retardos (tau) según min/max freq
    float maxF = YIN_MAX_FREQUENCY;
    if (maxF > 0.45f * rate) maxF = 0.45f * rate;     // guardia por Nyquist
    int tauMin = (int)floorf(rate / maxF);
    if (tauMin < 1) tauMin = 1;
    int tauMax = (int)ceilf(rate / YIN_MIN_FREQUENCY);

    // La ventana debe soportar el tau máximo
    if (tauMax + 1 >= N) {
        // Ventana demasiado corta para detectar hasta la frecuencia mínima
        return 0.0f;
    }

    // Buffers de trabajo
    float *diff = (float*)malloc((size_t)(tauMax + 1) * sizeof(float));
    float *cmnd = (float*)malloc((size_t)(tauMax + 1) * sizeof(float));
    if (!diff || !cmnd) { free(diff); free(cmnd); return 0.0f; }

    // 1) Difference function d(tau)
    diff[0] = 0.0f;
    for (int tau = 1; tau <= tauMax; ++tau) {
        double sum = 0.0;
        const int limit = N - tau;
        for (int j = 0; j < limit; ++j) {
            //const float delta = x[j] - x[j + tau]; // TODO acá cambiar el x por el circular buffer
            //const float delta = buffer[position - j] - buffer[position - j - tau]; // TODO acá cambiar el x por el circular buffer
            const float delta = buffer[position - N + j] - buffer[position - N + j + tau]; // TODO acá cambiar el x por el circular buffer
            sum += (double)delta * (double)delta;
        }
        diff[tau] = (float)sum;
    }

    // 2) CMND d'(tau) = d(tau) / ( (1/tau) * sum_{j=1..tau} d(j) )
    cmnd[0] = 1.0f;
    double acc = 0.0;
    for (int tau = 1; tau <= tauMax; ++tau) {
        acc += (double)diff[tau];
        cmnd[tau] = (acc == 0.0) ? 1.0f : (diff[tau] * (float)((double)tau / acc));
    }

    // 3) Buscar primer tau por debajo del umbral y avanzar al mínimo local
    const float thresh = YIN_THRESHOLD;
    int tau = -1;
    for (int t = tauMin; t <= tauMax; ++t) {
        if (cmnd[t] < thresh) {
            while (t + 1 <= tauMax && cmnd[t + 1] < cmnd[t]) ++t;
            tau = t;
            break;
        }
    }

    // Fallback: si nunca cruza el umbral, usar el mínimo global en el rango
    if (tau < 0) {
        float best = FLT_MAX;
        int bestT = -1;
        for (int t = tauMin; t <= tauMax; ++t) {
            if (cmnd[t] < best) { best = cmnd[t]; bestT = t; }
        }
        tau = bestT;
        if (tau < 0) { free(diff); free(cmnd); return 0.0f; }
        // Si la “confianza” es muy baja, puedes decidir devolver 0:
        // if (best > 0.9f) { free(diff); free(cmnd); return 0.0f; }
    }

    // 4) Interpolación parabólica alrededor del mínimo en CMND
    const float ym1 = (tau - 1 >= 1) ? cmnd[tau - 1] : cmnd[tau];
    const float y0  = cmnd[tau];
    const float yp1 = (tau + 1 <= tauMax) ? cmnd[tau + 1] : cmnd[tau];
    float offset = yin_parabolic_min_pos(ym1, y0, yp1);
    if (offset > 1.0f)  offset = 1.0f;
    if (offset < -1.0f) offset = -1.0f;

    const float tau_refined = (float)tau + offset;
    const float freq = (tau_refined > 0.0f) ? (rate / tau_refined) : 0.0f;

    free(diff);
    free(cmnd);
    return (freq > 0.0f) ? freq : 0.0f;
}