// assignment_6.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "stdlib.h"
#include "stdio.h"
#include "time.h"
#include "math.h"
#define P 12

long double universe_centroid[13];
long double vectors[6340][13];
long double code_book[8][13];
int code_book_size = 1;
long double distorsion_quatity = 0.03;
int number_of_vectors_in_universe = 6340;

long double W[13] = {0.0, 1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};

long double distorsion[8];
long double new_centroids[8][13];
int classified_vectors[8];

int read_vectors_into_array(int number_of_vectors, char *input_file)
{

    long double k;

    FILE *fp = fopen(input_file, "r");

    for (int i = 0; i < number_of_vectors; i++)
    {
        vectors[i][0] = 0;
        for (int j = 1; j <= 11; j++)
        {
            fscanf(fp, "%Lf,", &k);
            vectors[i][j] = k;
        }
        fscanf(fp, "%Lf\n", &k);
        vectors[i][12] = k;
    }

    return 0;
}

void find_centroid_of_universe()
{

    for (int i = 0; i < 13; i++)
        universe_centroid[i] = 0;

    for (int i = 0; i < 6340; i++)
    {
        for (int j = 1; j < 13; j++)
        {
            universe_centroid[j] += vectors[i][j];
        }
    }

    for (int i = 0; i < 13; i++)
    {
        universe_centroid[i] = (universe_centroid[i] / 6340);
    }

    for (int i = 0; i < 13; i++)
    {
        code_book[0][i] = universe_centroid[i];
    }
}

void split_vectors()
{
    long double temp_book[8][13];

    for (int i = 0; i < 8; i++)
    {
        for (int j = 0; j < 13; j++)
        {
            temp_book[i][j] = 0;
        }
    }

    int k = 0;
    for (int i = 0; i < code_book_size; i++)
    {
        for (int j = 1; j < 13; j++)
        {
            temp_book[k][j] = code_book[i][j] + distorsion_quatity;
            temp_book[k + 1][j] = code_book[i][j] - distorsion_quatity;
        }
        k = k + 2;
    }

    code_book_size = 2 * code_book_size;

    for (int i = 0; i < code_book_size; i++)
    {
        for (int j = 0; j < 13; j++)
        {
            code_book[i][j] = temp_book[i][j];
        }
    }
}

int find_minimum_index(long double *distance, int size)
{
    int minimum_index = 0;
    long double minimum = distance[0];
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

long double measure_tokuras_distance(long double Cr[13], long double Ct[13])
{
    long double distance = 0;
    for (int i = 1; i <= P; i++)
    {
        distance += W[i] * ((Cr[i] - Ct[i]) * (Cr[i] - Ct[i]));
    }
    return distance;
}

void classify_into_buckets(int number_of_vector_in_universe, int code_book_size)
{

    for (int i = 0; i < 8; i++)
    {
        distorsion[i] = 0;
        classified_vectors[i] = 0;
    }

    for (int i = 0; i < number_of_vector_in_universe; i++)
    {
        long double temp[8];

        for (int j = 0; j < code_book_size; j++)
        {

            temp[j] = measure_tokuras_distance(code_book[j], vectors[i]);
        }

        int bucket_index = find_minimum_index(temp, code_book_size);

        distorsion[bucket_index] += temp[bucket_index];

        new_centroids[bucket_index][0] = 0;

        for (int k = 1; k <= 12; k++)
        {
            new_centroids[bucket_index][k] += vectors[i][k];
        }

        classified_vectors[bucket_index]++;

        // printf("\nThe bucket index: %d\n", bucket_index);
    }

    for (int i = 0; i < code_book_size; i++)
    {
        for (int k = 1; k <= 12; k++)
        {
            new_centroids[i][k] = new_centroids[i][k] / classified_vectors[i];
        }
    }
}

int update_code_book(int code_book_size)
{
    for (int i = 0; i < code_book_size; i++)
    {
        for (int k = 1; k <= 12; k++)
        {
            code_book[i][k] = new_centroids[i][k];
        }
    }
    return 0;
}

long double get_distortion(int code_book_size)
{
    long double total_distortion = 0;
    for (int i = 0; i < code_book_size; i++)
    {
        total_distortion += distorsion[i];
    }
    return total_distortion;
}

void print_code_book(int code_book_size)
{
    printf("\n\nThe code book is :\n");
    for (int i = 0; i < code_book_size; i++)
    {
        for (int j = 0; j <= 12; j++)
        {
            printf("%Lf,", code_book[i][j]);
        }
        printf("\n");
    }
}

int _tmain(int argc, _TCHAR *argv[])
{
    long double total_distortion = 0, prev_distortion = 0;
    read_vectors_into_array(number_of_vectors_in_universe, "Universe.csv");
    find_centroid_of_universe();
    int m = 0;

    while (1)
    {
        split_vectors();

        do
        {
            prev_distortion = total_distortion;

            classify_into_buckets(number_of_vectors_in_universe, code_book_size);

            update_code_book(code_book_size);

            total_distortion = get_distortion(code_book_size);

            m++;
        } while (fabsl(prev_distortion - total_distortion) > 0.0001);

        if (code_book_size == 8)
        {
            break;
        }
    }
    print_code_book(code_book_size);

    // classify_into_buckets(number_of_vectors_in_universe, code_book_size);

    return 0;
}
