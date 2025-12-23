#include "math.h"
#include <sys/types.h>
#include <float.h>
#include <fftw3.h>

#include <fftw3.h>
#include <math.h>
#include <string.h>

#define F_MIN 80      // Hz (Mi2 ~82 Hz, nota más baja)
#define F_MAX 1200    // Hz (trastes altos, agudos)

// Función autocorrelacionada mejorada:
float fft_autocorr_freq(float *buffer, int position, int buffer_size, int rate) {
    int buffer_part = 16384; // o el que uses habitualmente

    // Asegúrate que puedes tomar buffer_part muestras...
    int start = (position > buffer_part) ? position - buffer_part : 0;

    // FFT directa
    fftwf_complex *out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * (buffer_part/2 + 1));
    fftwf_plan p = fftwf_plan_dft_r2c_1d(buffer_part, buffer + start, out, FFTW_ESTIMATE);
    fftwf_execute(p);

    // Espectro de potencia
    for (int i = 0; i < buffer_part/2 + 1; i++) {
        float re = out[i][0], im = out[i][1];
        out[i][0] = re*re + im*im;
        out[i][1] = 0.f;
    }

    // IFFT para autocorrelación
    float *autocorr = (float*) fftwf_malloc(sizeof(float) * buffer_part);
    fftwf_plan p2 = fftwf_plan_dft_c2r_1d(buffer_part, out, autocorr, FFTW_ESTIMATE);
    fftwf_execute(p2);

    // Normalización
    for (int i = 0; i < buffer_part; i++) autocorr[i] /= (float)buffer_part;

    // Búsqueda de lag útil (correspondiente a los periodos típicos)
    int lag_min = (int)(rate / F_MAX); // Lag máx frecuencia guitarra (~1200Hz)
    int lag_max = (int)(rate / F_MIN); // Lag mín frecuencia guitarra (~80Hz)
    if (lag_max > buffer_part-1) lag_max = buffer_part-1;

    // Salta a lag_min, busca máximo entre lag_min y lag_max
    float max_val = -1e10;
    int max_lag = lag_min;
    for (int i = lag_min; i <= lag_max; i++) {
        if (autocorr[i] > max_val) {
            max_val = autocorr[i];
            max_lag = i;
        }
    }

    // Liberar recursos
    fftwf_destroy_plan(p);
    fftwf_destroy_plan(p2);
    fftwf_free(out);
    fftwf_free(autocorr);

    // Evita división por cero
    if (max_lag == 0) return 0.0f;
    return (float)rate / (float)max_lag;
}

float fft_freq(float *buffer, int position, int buffer_size, int rate){
    // quizá se debiera reordenar el buffer...
    //int buffer_part = rate / 10;
    //int buffer_part =4096;  // 2^12
    //int buffer_part =8192;  // 2^13
    int buffer_part = 16384;  // 2^14
    fftwf_plan p;
    int create_wisdom = 1;
    int start;
    // Crear un array para almacenar el resultado de la FFT
    fftwf_complex *out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * (buffer_part + 1));

    if(position > buffer_part)
        start = position - buffer_part;
    else
        start = 0;
    
    // Crear un plan para FFTW
    if (fftwf_import_system_wisdom() != 0){
        // p = fftwf_plan_dft_r2c_1d(buffer_part, buffer + (buffer_part), out, FFTW_WISDOM_ONLY|FFTW_ESTIMATE);
        p = fftwf_plan_dft_r2c_1d(buffer_part, buffer + start, out, FFTW_WISDOM_ONLY|FFTW_ESTIMATE);
        printf("PitchDetection: using system wisdom file.\n");
    }
    else if (fftwf_import_wisdom_from_filename("wisdom") != 0)
	{
		//p = fftwf_plan_dft_r2c_1d(buffer_part, buffer + (buffer_part), out, FFTW_WISDOM_ONLY|FFTW_ESTIMATE);
        p = fftwf_plan_dft_r2c_1d(buffer_part, buffer + start, out, FFTW_WISDOM_ONLY|FFTW_ESTIMATE);
		//printf("PitchDetection: using plugin-provided wisdom file.\n");
	}
    else{
        if (create_wisdom){
            //fftwf_plan plan = fftwf_plan_dft_r2c_1d(buffer_part, buffer + (buffer_part), out, FFTW_MEASURE);
            fftwf_plan plan = fftwf_plan_dft_r2c_1d(buffer_part, buffer + start, out, FFTW_MEASURE);
            // Generar y guardar la wisdom
            if (fftwf_export_wisdom_to_filename("wisdom")) {
                printf("Wisdom file saved.\n");
                create_wisdom = 0;
            }
        }
        //p = fftwf_plan_dft_r2c_1d(buffer_part, buffer + (buffer_part), out, FFTW_ESTIMATE);
        p = fftwf_plan_dft_r2c_1d(buffer_part, buffer + start, out, FFTW_ESTIMATE);
        //printf("PitchDetection: failed to import wisdom file '%s', using estimate instead\n", wisdomFile);
        printf("PitchDetection: failed to import wisdom file, using estimate instead\n");
    }
    // Ejecutar la FFT
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