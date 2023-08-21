#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define P 12
#define FRAME_SIZE 320
#define FILES_NUM 10

long double *R;
long double *E_s;
long double *k_s;
long double *alpha_s;
long double *A_s;
long double *Ci_s;

int speech_data[100000];
int amplitude_count = 0;

long double k(int i);
long double E(int i);
long double alpha(int x, int p);

/*
    Multiplying the amplitude with hamming window
*/

void multiply_with_hamming_window(long double *amplitudes, int size)
{
    for (int i = 0; i < size; i++)
    {
        long double hamming_window = 0.54 - 0.46 * (((2 * 3.141592 * i) / (size - 1)));
        printf("%Lf\n", hamming_window);
        amplitudes[i] = (amplitudes[i] * hamming_window);
    }
}

/*
    Calcuting Ri's part of durbin's algorithm
*/

void calculate_Ris(long double *amplitude_array, int size)
{
    R = (long double *)calloc((P + 1), sizeof(long double));
    multiply_with_hamming_window(amplitude_array, size);

    long double sum_of_shifted_amplitudes = 0;
    for (int shift = 0; shift <= P; shift++)
    {
        for (int i = 0; i < size - shift; i++)
        {
            long double squared_value = amplitude_array[i] * amplitude_array[i + shift];

            sum_of_shifted_amplitudes += squared_value;
        }
        R[shift] = sum_of_shifted_amplitudes / (long double)FRAME_SIZE;
        sum_of_shifted_amplitudes = 0;
    }
}

/*
    For making the amplitudes normalized
*/

long double *normalize_data(long double *frame_amplitudes, long double max, long double scale)
{

    long double *normalize_frames = (long double *)malloc(FRAME_SIZE * sizeof(long double));
    for (int i = 0; i < FRAME_SIZE; i++)
    {
        long double a = ((frame_amplitudes[i] / abs(max)) * scale);
        normalize_frames[i] = a;
    }

    return normalize_frames;
}

long double E(int i)
{
    if (E_s[i] != 0)
        return E_s[i];

    if (i == 0)
        return R[0];
    long double temp = k(i);

    temp = temp * temp;

    long double temp1 = E(i - 1);
    E_s[i] = ((1 - temp) * temp1);

    return E_s[i];
}

long double k(int i)
{
    if (k_s[i] != 0)
        return k_s[i];

    long double temp = R[i];
    long double sum = 0;

    for (int j = 1; j <= i - 1; j++)
    {
        long double temp1 = alpha(j, i - 1);
        long double temp2 = R[i - j];
        sum += (temp1 * temp2);
    }

    long double temp3 = E(i - 1);
    long double result = (temp - sum) / temp3;

    k_s[i] = result;

    return result;
}

long double alpha(int x, int p)
{
    if (x == p)
        return k(x);
    else
    {
        long double temp = alpha(x, p - 1);
        long double temp1 = k(p);
        long double temp2 = alpha(p - x, p - 1);
        return temp - temp1 * temp2;
    }
}

/*
    Calcuting Ai's part of durbin's algorithm
*/

void calculate_Ais()
{
    A_s = (long double *)calloc(P + 1, sizeof(long double));
    E_s = (long double *)calloc((P + 1), sizeof(long double));
    k_s = (long double *)calloc(P + 1, sizeof(long double));
    for (int i = 1; i <= P; i++)
    {
        A_s[i] = alpha(i, P);
    }
}

void calculate_Cis()
{
    Ci_s = (long double *)calloc(P + 1, sizeof(long double));
    Ci_s[0] = (R[0]);

    for (int m = 1; m <= P; m++)
    {
        long double temp_sum = 0;

        for (int k = 1; k < m; k++)
        {
            temp_sum += ((long double)k / (long double)m) * Ci_s[k] * A_s[m - k];
        }

        Ci_s[m] = A_s[m] + temp_sum;
    }
}

long double read_data_and_find_max(char *input_file)
{
    FILE *ptr = fopen(input_file, "r");
    long double max = 0;
    char line[10];
    amplitude_count = 0;

    if (NULL == ptr)
    {
        printf("file can't be opened \n");
        return 0;
    }

    while (fgets(line, sizeof(line), ptr) != NULL)
    {
        long double amplitude = atof(line);
        speech_data[amplitude_count++] = amplitude;
        if (fabs(amplitude) > max)
        {
            max = fabs(amplitude);
        }
    }
    fclose(ptr);

    return max;
}

long double *get_amplitudes_in_frame(int start, int end)
{
    int count = 0;
    long double *frame = (long double *)malloc(320 * sizeof(long double));
    for (int i = start; i < end; i++)
    {
        frame[count++] = speech_data[i];
    }
    return frame;
}

void write_to_universe_file(long double *Ci_s)
{
    FILE *ptr = fopen("Universe.txt", "a+");
    for (int i = 1; i < P; i++)
    {
        fprintf(ptr, "%Lf,", Ci_s[i]);
    }
    fprintf(ptr, "%Lf\n", Ci_s[P]);
    fclose(ptr);
}

void read_data_and_generate_cepstrals(int amplitude_count, long double max)
{
    int start = 0;
    int end = 320;

    while (end < amplitude_count)
    {
        long double *frame = get_amplitudes_in_frame(start, end);
        frame = normalize_data(frame, max, 5000);
        calculate_Ris(frame, FRAME_SIZE);
        calculate_Ais();
        calculate_Cis();
        write_to_universe_file(Ci_s);
        free(R);
        free(E_s);
        free(k_s);
        free(alpha_s);
        free(A_s);
        free(Ci_s);
        start = start + 80;
        end = end + 80;
    }
}

int main()
{

    long double max = read_data_and_find_max("1/224101014_E_1_5.txt");
    read_data_and_generate_cepstrals(amplitude_count, max);

    return 0;
}