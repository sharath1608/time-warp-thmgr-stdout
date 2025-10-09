// Dynamic Time Warping (DTW) algorithm, which is widely used in speech recognition, time series analysis, and pattern matching.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>

#define MAX_SEQUENCE_LENGTH 10000
#define WINDOW_SIZE 100  // Sakoe-Chiba band width

// Structure to store DTW result
typedef struct {
    double distance;
    int* path_i;
    int* path_j;
    int path_length;
} DTWResult;

// Function to generate random sequence
void generate_sequence(double* sequence, int length) {
    for(int i = 0; i < length; i++) {
        sequence[i] = ((double)rand() / RAND_MAX) * 100.0;
    }
}

// Function to calculate Euclidean distance
double calculate_distance(double a, double b) {
    return (a - b) * (a - b);
}

// Function to find minimum of three values
double min3(double a, double b, double c) {
    return fmin(fmin(a, b), c);
}

// Main DTW implementation with windowing
DTWResult* compute_dtw(double* sequence1, int len1, 
                      double* sequence2, int len2) {
    // Allocate cost matrix
    double** cost_matrix = (double**)malloc(len1 * sizeof(double*));
    for(int i = 0; i < len1; i++) {
        cost_matrix[i] = (double*)malloc(len2 * sizeof(double));
        for(int j = 0; j < len2; j++) {
            cost_matrix[i][j] = DBL_MAX;
        }
    }

    // Initialize first element
    cost_matrix[0][0] = calculate_distance(sequence1[0], sequence2[0]);

    // Initialize first row and column
    for(int i = 1; i < len1; i++) {
        if(i <= WINDOW_SIZE) {
            cost_matrix[i][0] = cost_matrix[i-1][0] + 
                               calculate_distance(sequence1[i], sequence2[0]);
        }
    }
    for(int j = 1; j < len2; j++) {
        if(j <= WINDOW_SIZE) {
            cost_matrix[0][j] = cost_matrix[0][j-1] + 
                               calculate_distance(sequence1[0], sequence2[j]);
        }
    }

    // Fill cost matrix with windowing (Sakoe-Chiba band)
    for(int i = 1; i < len1; i++) {
        int start = (i - WINDOW_SIZE > 1) ? i - WINDOW_SIZE : 1;
        int end = (i + WINDOW_SIZE < len2) ? i + WINDOW_SIZE : len2;
        
        for(int j = start; j < end; j++) {
            double cost = calculate_distance(sequence1[i], sequence2[j]);
            cost_matrix[i][j] = cost + min3(
                cost_matrix[i-1][j],    // insertion
                cost_matrix[i][j-1],    // deletion
                cost_matrix[i-1][j-1]   // match
            );
        }
    }

    // Backtracking to find optimal path
    DTWResult* result = (DTWResult*)malloc(sizeof(DTWResult));
    result->path_i = (int*)malloc((len1 + len2) * sizeof(int));
    result->path_j = (int*)malloc((len1 + len2) * sizeof(int));
    
    int path_index = 0;
    int i = len1 - 1;
    int j = len2 - 1;
    
    while(i > 0 || j > 0) {
        result->path_i[path_index] = i;
        result->path_j[path_index] = j;
        path_index++;
        
        if(i == 0) {
            j--;
        } else if(j == 0) {
            i--;
        } else {
            double min_cost = min3(
                cost_matrix[i-1][j],
                cost_matrix[i][j-1],
                cost_matrix[i-1][j-1]
            );
            
            if(min_cost == cost_matrix[i-1][j-1]) {
                i--; j--;
            } else if(min_cost == cost_matrix[i-1][j]) {
                i--;
            } else {
                j--;
            }
        }
    }
    
    // Add final point
    result->path_i[path_index] = 0;
    result->path_j[path_index] = 0;
    path_index++;
    
    result->path_length = path_index;
    result->distance = cost_matrix[len1-1][len2-1];

    // Free cost matrix
    for(int i = 0; i < len1; i++) {
        free(cost_matrix[i]);
    }
    free(cost_matrix);

    return result;
}

// Function to perform subsequence DTW search
int find_best_match(double* long_sequence, int long_length,
                   double* pattern, int pattern_length,
                   double* best_distance) {
    int best_position = 0;
    *best_distance = DBL_MAX;

    // Sliding window approach
    for(int i = 0; i <= long_length - pattern_length; i++) {
        // Extract subsequence
        double* subsequence = (double*)malloc(pattern_length * sizeof(double));
        for(int j = 0; j < pattern_length; j++) {
            subsequence[j] = long_sequence[i + j];
        }

        // Compute DTW
        DTWResult* result = compute_dtw(subsequence, pattern_length,
                                      pattern, pattern_length);

        // Update best match if necessary
        if(result->distance < *best_distance) {
            *best_distance = result->distance;
            best_position = i;
        }

        // Cleanup
        free(subsequence);
        free(result->path_i);
        free(result->path_j);
        free(result);
    }

    return best_position;
}

int main(int argc, char * argv[]) {

    srand(time(NULL));

    // Generate sample sequences
    int long_length = atoi(argv[1]);
    int pattern_length = 500;
    
    double* long_sequence = (double*)malloc(long_length * sizeof(double));
    double* pattern = (double*)malloc(pattern_length * sizeof(double));
    
    generate_sequence(long_sequence, long_length);
    generate_sequence(pattern, pattern_length);

    // Measure execution time
    clock_t start = clock();
    
    // Find best matching subsequence
    double best_distance;
    int best_position = find_best_match(long_sequence, long_length,
                                      pattern, pattern_length,
                                      &best_distance);
    
    clock_t end = clock();
    double cpu_time = ((double) (end - start)) / CLOCKS_PER_SEC;

    // Print results
    printf("Best match found at position: %d\n", best_position);
    printf("Distance: %f\n", best_distance);
    printf("Execution time: %f seconds\n", cpu_time);

    // Cleanup
    free(long_sequence);
    free(pattern);

    return 0;
}
