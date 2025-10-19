#include "math.h"
#include <sys/types.h>
#include <float.h>

#include <fftw3.h>

//ajustar variables de input
float corr_freq(float *buffer, int position, int buffer_size, int rate){
//float corr_freq(float* buffer, int position, float sample_rate) {
    float max_corr = 0.0f;
    int max_index = 0;

    for (int delay = 1; delay < buffer_size; delay++) {
        float corr = 0.0f;
        for (int i = 0; i < position - delay; i++) {
            corr += buffer[i] * buffer[i + delay];
        }

        // Normaliza la correlación
        corr /= (position - delay);

        if (corr > max_corr) {
            max_corr = corr;
            max_index = delay;
        }
    }

    // Calcular la frecuencia fundamental
    return rate / max_index;
}


float fft_freq(float *buffer, int position, int buffer_size, int rate){
    // quizá se debiera reordenar el buffer...
    // Crear el plan FFTW
    int buffer_part = buffer_size / 10;
    fftwf_complex *out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * (buffer_part + 1));
    fftwf_plan p = fftwf_plan_dft_r2c_1d(buffer_part, buffer + (buffer_part), out, FFTW_ESTIMATE);
    
    // Realizar la FFT
    fftwf_execute(p);   
    // Encontrar la frecuencia máxima
    double max_amplitude = 0;
    int max_index = 0;
    
    for (int i = 0; i < buffer_part + 1; i++) {
        double amplitude = sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]);
        if (amplitude > max_amplitude) {
            max_amplitude = amplitude;
            max_index = i;
        }
    }
    
    // Liberar recursos
    fftwf_destroy_plan(p);
    fftwf_free(out);
    return(max_index * rate / buffer_part);
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


#define MAX_LAG 300
#define THRESHOLD 0.1
float yin_freq(float* buffer, int position, int buffer_size, int rate){
    if (!buffer || buffer_size < 4 || rate <= 0) return 0.0f;
    float difference[MAX_LAG] = {0};
    float cumulative[MAX_LAG] = {0};
    int tau, i;
    
    // Paso 1: Función de diferencia
    for (tau = 1; tau < MAX_LAG; tau++) {
        for (i = 0; i < buffer_size - tau; i++) {
            float diff = buffer[i] - buffer[i + tau];
            difference[tau] += diff * diff;
        }
        difference[tau] /= (buffer_size - tau);
    }
    
    // Paso 2: Función acumulativa
    cumulative[0] = 1.0;
    float sum = 0.0;
    for (tau = 1; tau < MAX_LAG; tau++) {
        sum += difference[tau];
        cumulative[tau] = difference[tau] / (sum / tau);
    }
    
    // Paso 3: Buscar el mínimo por debajo del umbral
    int min_tau = -1;
    for (tau = 2; tau < MAX_LAG; tau++) {
        if (cumulative[tau] < THRESHOLD) {
            min_tau = tau;
            break;
        }
    }
    
    // Paso 4: Convertir a frecuencia
    if (min_tau != -1) {
        return rate / min_tau;
    } else {
        return 0.0; // No se detectó pitch
    } 
}