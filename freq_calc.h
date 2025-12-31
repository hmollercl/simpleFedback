#include "math.h"
#include <sys/types.h>
#include <float.h>
#include <fftw3.h>

#define F_MIN 300.0f     // Hz (Mi2 ~82 Hz, nota más baja)
#define F_MAX 1200.0f    // Hz (trastes altos, agudos)
#define FFT_SIZE 16384   // tamaño de la ventana para calcular fft



// RT-safe: NO malloc, NO printf, NO creación/destrucción de planes.
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
    for (int i = 0; i < N; ++i) {
        m->pd_time[i] = buffer[idx];
        idx++;
        if (idx >= buffer_size) {
            idx = 0;
        }
    }

    // === 2. FFT (forward) ===
    fftwf_execute(m->pd_plan_fwd);

    // === 3. Espectro de potencia: |X(f)|^2 ===
    int spec_size = N / 2 + 1;
    for (int i = 0; i < spec_size; ++i) {
        float re = m->pd_freq[i][0];
        float im = m->pd_freq[i][1];
        m->pd_freq[i][0] = re * re + im * im; // potencia real
        m->pd_freq[i][1] = 0.0f;              // parte imaginaria = 0
    }

    // === 4. IFFT → autocorrelación ===
    fftwf_execute(m->pd_plan_inv);

    // === 5. Normalizar autocorrelación ===
    for (int i = 0; i < N; ++i) {
        m->pd_autocorr[i] /= (float)N;
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

    // === 7. Convertir lag a frecuencia ===
    float freq = (float)m->rate / (float)max_lag;

    // (Opcional) puedes hacer un pequeño suavizado aquí usando m->calc_freq previo:
    // float alpha = 0.3f;
    // freq = alpha * freq + (1.0f - alpha) * m->calc_freq;

    return freq;
}


// Función autocorrelacionada mejorada:
float fft_autocorr_freq(float *buffer, int position, int buffer_size, int rate) {
    int buffer_part = 16384; // o el que uses habitualmente
    int start;
    fftwf_plan p;
    fftwf_plan p2;

    // Asegúrate que puedes tomar buffer_part muestras...
    start = (position > buffer_part) ? position - buffer_part : 0;

    // FFT directa
    fftwf_complex *out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * (buffer_part/2 + 1));
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
        p = fftwf_plan_dft_r2c_1d(buffer_part, buffer + start, out, FFTW_ESTIMATE);
        printf("PitchDetection: using estimate.\n");
    }
    
    fftwf_execute(p);

    // Espectro de potencia
    for (int i = 0; i < buffer_part/2 + 1; i++) {
        float re = out[i][0], im = out[i][1];
        out[i][0] = re*re + im*im;
        out[i][1] = 0.f;
    }

    // IFFT para autocorrelación
    float *autocorr = (float*) fftwf_malloc(sizeof(float) * buffer_part);
    /*if (fftwf_import_system_wisdom() != 0){
        p2 = fftwf_plan_dft_r2c_1d(buffer_part, out, autocorr, FFTW_WISDOM_ONLY|FFTW_ESTIMATE);
        printf("PitchDetection: using system wisdom file.\n");
    }
    else if (fftwf_import_wisdom_from_filename("wisdom") != 0)
	{

        p2 = fftwf_plan_dft_r2c_1d(buffer_part, out, autocorr, FFTW_WISDOM_ONLY|FFTW_ESTIMATE);
		//printf("PitchDetection: using plugin-provided wisdom file.\n");
	}
    else{
        p2 = fftwf_plan_dft_c2r_1d(buffer_part, out, autocorr, FFTW_ESTIMATE);
        printf("PitchDetection: using estimate.\n");
    }*/
    p2 = fftwf_plan_dft_c2r_1d(buffer_part, out, autocorr, FFTW_ESTIMATE);
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
    //int buffer_part =4096;  // 2^12
    //int buffer_part =8192;  // 2^13
    int buffer_part = 16384;  // 2^14
    fftwf_plan p;
    int create_wisdom = 1;
    int start;
    // Crear un array para almacenar el resultado de la FFT
    fftwf_complex *out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * (buffer_part/2 + 1));

    /*if(position > buffer_part)
        start = position - buffer_part;
    else
        start = 0;*/
    start = (position > buffer_part) ? position - buffer_part : 0;
    
    // Crear un plan para FFTW
    if (fftwf_import_system_wisdom() != 0){
        p = fftwf_plan_dft_r2c_1d(buffer_part, buffer + start, out, FFTW_WISDOM_ONLY|FFTW_ESTIMATE);
        printf("PitchDetection: using system wisdom file.\n");
    }
    else if (fftwf_import_wisdom_from_filename("wisdom") != 0)
	{
        p = fftwf_plan_dft_r2c_1d(buffer_part, buffer + start, out, FFTW_WISDOM_ONLY|FFTW_ESTIMATE);
		//printf("PitchDetection: using plugin-provided wisdom file.\n");
	}
    else{
        if (create_wisdom){
            fftwf_plan plan = fftwf_plan_dft_r2c_1d(buffer_part, buffer + start, out, FFTW_MEASURE);
            // Generar y guardar la wisdom
            if (fftwf_export_wisdom_to_filename("wisdom")) {
                printf("Wisdom file saved.\n");
                create_wisdom = 0;
            }
        }
        p = fftwf_plan_dft_r2c_1d(buffer_part, buffer + start, out, FFTW_ESTIMATE);
        //printf("PitchDetection: failed to import wisdom file '%s', using estimate instead\n", wisdomFile);
        printf("PitchDetection: failed to import wisdom file, using estimate instead\n");
    }
    // Ejecutar la FFT
    fftwf_execute(p);   

    // Encontrar la frecuencia máxima
    double max_amplitude = 0;
    int max_index = 0;
    
    for (int i = 0; i < buffer_part/2 + 1; i++) {
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