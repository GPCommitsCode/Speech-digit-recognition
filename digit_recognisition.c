// assignment_4.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#define P 12
#define FRAME_SIZE 320
#define FILES_NUM 10

double long *R;
double long *E_s;
double long *k_s;
double long *alpha_s;
double long *A_s;
double long *Ci_s;
double long Energy[100];

// We are taking 13 values for cepstral array for starting index from 1

double long Ci_s_list[100][13]; // for storing the cepstal coefficients for one file
double long steady_Cis[5][13];  // for storing the steady cepstral coefficients for one file
double long populated_code_book[5][13];

double long k(int i);
double long E(int i);
int frames_in_current_file = 0;

// Vowels listed in the order for classification ease
char vowels[5] = {'a', 'e', 'i', 'o', 'u'};

// The reference file names
char file_name[5][10][15] = {
    {"a_1_ref.txt", "a_2_ref.txt", "a_3_ref.txt", "a_4_ref.txt", "a_5_ref.txt", "a_6_ref.txt", "a_7_ref.txt", "a_8_ref.txt", "a_9_ref.txt", "a_10_ref.txt"},
    {"e_1_ref.txt", "e_2_ref.txt", "e_3_ref.txt", "e_4_ref.txt", "e_5_ref.txt", "e_6_ref.txt", "e_7_ref.txt", "e_8_ref.txt", "e_9_ref.txt", "e_10_ref.txt"},
    {"i_1_ref.txt", "i_2_ref.txt", "i_3_ref.txt", "i_4_ref.txt", "i_5_ref.txt", "i_6_ref.txt", "i_7_ref.txt", "i_8_ref.txt", "i_9_ref.txt", "i_10_ref.txt"},
    {"o_1_ref.txt", "o_2_ref.txt", "o_3_ref.txt", "o_4_ref.txt", "o_5_ref.txt", "o_6_ref.txt", "o_7_ref.txt", "o_8_ref.txt", "o_9_ref.txt", "o_10_ref.txt"},
    {"u_1_ref.txt", "u_2_ref.txt", "u_3_ref.txt", "u_4_ref.txt", "u_5_ref.txt", "u_6_ref.txt", "u_7_ref.txt", "u_8_ref.txt", "u_9_ref.txt", "u_10_ref.txt"}};

// The test file names

char test_file_name[5][10][15] = {
    {"a_1_test.txt", "a_2_test.txt", "a_3_test.txt", "a_4_test.txt", "a_5_test.txt", "a_6_test.txt", "a_7_test.txt", "a_8_test.txt", "a_9_test.txt", "a_10_test.txt"},
    {"e_1_test.txt", "e_2_test.txt", "e_3_test.txt", "e_4_test.txt", "e_5_test.txt", "e_6_test.txt", "e_7_test.txt", "e_8_test.txt", "e_9_test.txt", "e_10_test.txt"},
    {"i_1_test.txt", "i_2_test.txt", "i_3_test.txt", "i_4_test.txt", "i_5_test.txt", "i_6_test.txt", "i_7_test.txt", "i_8_test.txt", "i_9_test.txt", "i_10_test.txt"},
    {"o_1_test.txt", "o_2_test.txt", "o_3_test.txt", "o_4_test.txt", "o_5_test.txt", "o_6_test.txt", "o_7_test.txt", "o_8_test.txt", "o_9_test.txt", "o_10_test.txt"},
    {"u_1_test.txt", "u_2_test.txt", "u_3_test.txt", "u_4_test.txt", "u_5_test.txt", "u_6_test.txt", "u_7_test.txt", "u_8_test.txt", "u_9_test.txt", "u_10_test.txt"}};

// Tokura's weight given by sir, after years of experimentation

double long W[13] = {0.0, 1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};

double long alpha(int x, int p)
{
        if (x == p)
        return k(x);
        else    
    {
                double long temp = alpha(x, p - 1);
                double long temp1 = k(p);
                double long temp2 = alpha(p - x, p - 1);
                return temp - temp1 * temp2;
           
    }
}

double long k(int i)
{
        if (k_s[i] != 0)
        return k_s[i];

        double long temp = R[i];
        double long sum = 0;

        for (int j = 1; j <= i - 1; j++)
   
    {
                double long temp1 = alpha(j, i - 1);
                double long temp2 = R[i - j];
                sum += (temp1 * temp2);
           
    }

        double long temp3 = E(i - 1);
        double long result = (temp - sum) / temp3;

        k_s[i] = result;

        return result;
}

double long E(int i)
{
        if (E_s[i] != 0)
        return E_s[i];

        if (i == 0)
        return R[0];
        double long temp = k(i);

        temp = temp * temp;

        double long temp1 = E(i - 1);
        E_s[i] = ((1 - temp) * temp1);

        return E_s[i];
}

/*
    Multiplying the amplitude with hamming window
*/

void multiply_with_hamming_window(long double *amplitudes, int size)
{
        // printf("\nhamming window:\n");
    for (int i = 0; i < size; i++)
   
    {
                long double hamming_window = 0.54 - 0.46 * (cosl((2 * 3.141592 * i) / (size - 1)));
              //  printf("%Lf\n", hamming_window);
        amplitudes[i] = (amplitudes[i] * hamming_window);
           
    }
}

/*
    Calcuting Ri's part of durbin's algorithm
*/

void calculate_Ris(double long *amplitude_array, int size)
{
        R = (double long *)calloc((P + 1), sizeof(double long));
        multiply_with_hamming_window(amplitude_array, size);

        long double sum_of_shifted_amplitudes = 0;
        for (int shift = 0; shift <= P; shift++)
   
    {
                for (int i = 0; i < size - shift; i++)
       
        {
                        long double squared_value = amplitude_array[i] * amplitude_array[i + shift];

                        sum_of_shifted_amplitudes += squared_value;
                   
        }
                R[shift] = sum_of_shifted_amplitudes / (double long)FRAME_SIZE;
                sum_of_shifted_amplitudes = 0;
           
    }
}

/*
    Calcuting Ai's part of durbin's algorithm
*/

void calculate_Ais()
{
        A_s = (double long *)calloc(P + 1, sizeof(double long));
        E_s = (double long *)calloc((P + 1), sizeof(double long));
        k_s = (double long *)calloc(P + 1, sizeof(double long));
        for (int i = 1; i <= P; i++)
   
    {
                A_s[i] = alpha(i, P);
           
    }
}

/*
    For finding the Ci's
*/

void calculate_Cis()
{
        Ci_s = (double long *)calloc(P + 1, sizeof(double long));
        Ci_s[0] = log(R[0]);

        for (int m = 1; m <= P; m++)
   
    {
                double long temp_sum = 0;

                for (int k = 1; k < m; k++)
       
        {
                        temp_sum += ((double long)k / (double long)m) * Ci_s[k] * A_s[m - k];
                   
        }

                Ci_s[m] = A_s[m] + temp_sum;
           
    }

      //  for (int i = 1; i <= P; i++)
  //  {
  //      printf("\nC%d: %Lf\n", i, Ci_s[i]);
  //  }
}

/*
    Calculate STE
*/

double long calculate_STE(double long *amplitude_array, int size)
{
    double long SUM_STE = 0;

    for (int i = 0; i < size; i++)
    {
        SUM_STE += amplitude_array[i] * amplitude_array[i];
    }

    return (SUM_STE / size);
}

/*
Finding the max from the input file

The value max is required for normalization
*/
double long find_max(char *input_file)
{
    FILE *ptr = fopen(input_file, "r");
    double long max = 0;
    char line[10];

    if (NULL == ptr)
    {
            printf("file can't be opened \n");
        return 0;
           
    }

    while (fgets(line, sizeof(line), ptr) != NULL)
    {
        double long amplitude = atof(line);
        if (fabs(amplitude) > max)
        {
            max = fabs(amplitude);
        }
    }
    fclose(ptr);

    return max;
}

/*
    For making the amplitudes normalized
*/

double long *normalize_data(double long *frame_amplitudes, double long max, double long scale)
{

    double long *normalize_frames = (double long *)malloc(FRAME_SIZE * sizeof(double long));
    for (int i = 0; i < FRAME_SIZE; i++)
    {
        double long a = ((frame_amplitudes[i] / abs(max)) * scale);
        normalize_frames[i] = a;
    }

    return normalize_frames;
}

/*
    Read data from file and transform them into frame, calculate also the energy for comparision
*/

void readDataFromFile(char *file_path, double long max, double long scale, int *frame_in_files)
{

    FILE *ptr = fopen(file_path, "r");
    double long amplitude_array[FRAME_SIZE];
    double long *normalized_amplitude_array;
    char line[10];
    int amp_count = 0;

    if (NULL == ptr)
    {
        printf("file can't be opened \n");
        return;
    }

    while (fgets(line, sizeof(line), ptr) != NULL)
    {
        amplitude_array[amp_count++] = atof(line);
        if (amp_count == FRAME_SIZE)
        {
            normalized_amplitude_array = normalize_data(amplitude_array, max, scale);
            calculate_Ris(normalized_amplitude_array, FRAME_SIZE);
            calculate_Ais();
            calculate_Cis();
            for (int i = 1; i <= P; i++)
            {
                Ci_s_list[(*frame_in_files)][i] = Ci_s[i]; // part where we calculate the cepstral coeffecients
            }
            double long energy = calculate_STE(normalized_amplitude_array, FRAME_SIZE);
            Energy[*frame_in_files] = energy; // calculate the corresponding energy

            (*frame_in_files)++;
            amp_count = 0;
        }
    }
}

/*
    For collecting steady part of the file - 5 frames - 2 left of the highesh energy frame + 2 right of the highest energy frame + 1 highest energy frame
*/

void get_steady_part(int frames_in_current_file, double long *Energy)
{
    int max_index = 0;
    int max = Energy[0];
    for (int i = 1; i < frames_in_current_file; i++)
    {
        if (max < Energy[i])
        {
            max = Energy[i];
            max_index = i;
        }
    }

    if (max_index - 2 >= 0 && max_index + 2 < FRAME_SIZE)
    {
        int count = 0;
        for (int i = max_index - 2; i <= max_index + 2; i++)
        {
            for (int k = 1; k <= P; k++)
            {
                steady_Cis[count][k] = Ci_s_list[i][k];
            }
            count++;
        }
    }
    else
    {
        printf("problem finding steady part!!");
    }
}

/*
    To measure the tokura's distance
*/
double long measure_tokuras_distance(double long Cr[13], double long Ct[13])
{
    double long distance = 0;
    for (int i = 1; i <= P; i++)
    {
        distance += W[i] * ((Cr[i] - Ct[i]) * (Cr[i] - Ct[i]));
    }
    return distance;
}

/*
    To populate code book from the file
*/
void populate_code_book(char *code_book_file_path)
{
    FILE *ptr = fopen(code_book_file_path, "r");
    double long k = 0;
    for (int i = 0; i < 5; i++)
    {
        populated_code_book[i][0] = 0;
        for (int j = 1; j <= 11; j++)
        {
            fscanf(ptr, "%Lf,", &k);
            populated_code_book[i][j] = k;
        }
        fscanf(ptr, "%Lf\n", &k);
        populated_code_book[i][12] = k;
    }
}

/*
    Finding the index which has minimum distance
*/
int find_minimum_index(double long *distance, int size)
{
    int minimum_index = 0;
    double long minimum = distance[0];
    for (int i = 1; i < size; i++)
    {
        if (distance[i] < minimum)
        {
            minimum = distance[i];
            minimum_index = i;
        }
    }
    return minimum_index;
}

/*
    Finding the highest frequency index
*/
int find_highest_frequency_index(int *frequency_array, int size)
{
    int highest_frequency_index = 0;
    int highest_frequency = 0;

    for (int i = 1; i < size; i++)
    {
        if (frequency_array[i] > highest_frequency)
        {
            highest_frequency = frequency_array[i];
            highest_frequency_index = i;
        }
    }
    return highest_frequency_index;
}

int main()
{

    int accurate_vowel_count = 0;

    // Code book : ---------------- This part can be run to generate code book--------------------------------
    /*
    for(int q = 0; q < 5; q++){

        double long *final_cepstral = (double long *)calloc(P+1, sizeof(double long));

        for(int i =0; i < FILES_NUM; i++){
            char *filename = file_name[q][i];

            double long max_in_file = find_max(filename);
            readDataFromFile(filename, max_in_file, 5000, &frames_in_current_file);
            get_steady_part(frames_in_current_file, Energy);
            for(int k = 0; k < 5; k++){
                for(int i = 1; i <=P; i++){
                    final_cepstral[i] += steady_Cis[k][i];
                }
            }

            frames_in_current_file = 0;
        }

        FILE *codebook = fopen("code_book.txt", "a");

        for(int i = 1; i <P; i++){
            fprintf(codebook, "%Lf,", (final_cepstral[i]/50));
        }

        fprintf(codebook, "%Lf\n",  (final_cepstral[P]/50));

        fclose(codebook);

    }
*/

    //--------------------------------------End of code book generation part--------------------------------------

    populate_code_book("code_book.txt"); // For populating the code book

    for (int j = 0; j < 5; j++)
    {
        for (int i = 0; i < FILES_NUM; i++)
        {

            //------------------ Part for reading from the file and getting the steady part(5 frames into) ----------------------------------
            char *filename = test_file_name[j][i];
            printf("\n%s:\n", filename);

            double long max_in_file = find_max(filename);
            readDataFromFile(filename, max_in_file, 5000, &frames_in_current_file);
            get_steady_part(frames_in_current_file, Energy);
            // ----------------- Ending of the common area for all the methods ---------------------------------------

            //--------------Average tokura's distance method : Take average take tokura's distance then classify----------------------------

            double long *final_cepstral = (double long *)calloc(P + 1, sizeof(double long));

            final_cepstral[0] = 0;

            for (int k = 0; k < 5; k++)
            {
                for (int i = 1; i <= P; i++)
                {
                    final_cepstral[i] += steady_Cis[k][i];
                }
            }

            // printf("\nFinal cepstrum is:\n");

            for (int i = 1; i <= P; i++)
            {
                final_cepstral[i] = (final_cepstral[i] / 5);
                // printf("%Lf\t", final_cepstral[i]);
            }

            // printf("\n");

            double long tokura_distance[5];

            for (int i = 0; i < 5; i++)
            {
                tokura_distance[i] = measure_tokuras_distance(populated_code_book[i], final_cepstral);
                // printf("\n%c: %Lf\n",vowels[i],tokura_distance[i]);
            }

            int vowel_index = find_minimum_index(tokura_distance, 5);
            printf("The vowel is %c\n", vowels[vowel_index]);
            if (vowel_index == j)
                accurate_vowel_count++;

            // ------------------------------------End of Average tokura's distance method---------------------------------------

            //----------------------------------------Frequency counting method--------------------------------------------------
            /*
                printf("\n---------------\n");

                double long tokura_distance[5];
                int *vowel_frequency_count = (int *)calloc(5, sizeof(int));
                for(int j=0; j<5; j++) {
                    for(int i=0; i<5; i++){
                        tokura_distance[i] = measure_tokuras_distance(populated_code_book[i],steady_Cis[j]);
                        //printf("\n%c: %Lf\n",vowels[i],tokura_distance[i]);
                    }

                    int vowel_index = find_minimum_index(tokura_distance, 5);
                    printf("The vowel is %c\n", vowels[vowel_index]);
                    vowel_frequency_count[vowel_index]++;
                }

                int classification_vowel_index = find_highest_frequency_index(vowel_frequency_count,5);
                printf("The highest frequncy is for %c so the file can be classified as %c file\n", vowels[classification_vowel_index], vowels[classification_vowel_index]);
                if(classification_vowel_index == j) accurate_vowel_count++;
                */

            //-----------------------------------------End of frequency counting method-----------------------------

            //-----------------------------------------start of ending common area -------------------

            frames_in_current_file = 0;

            printf("\n\n");
            // End
        }
    }

    printf("\nAccuracy: (correct_prediction_of_file)/(total_files(50)) is : %f percentage \n", ((float)accurate_vowel_count / (float)50) * 100);

        return 0;
}