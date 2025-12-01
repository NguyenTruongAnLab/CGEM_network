#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <nlopt.h>
#include <string.h>
#include <time.h>

#define NUM_PARAMS 6
#define NUM_DAYS 730
#define NUM_CELLS 200

int iteration_count = 0;
const int max_iterations = 10000; // Higher iteration count for complexity
double observed_data[NUM_DAYS][NUM_CELLS + 1]; // +1 for tidal range
FILE *output_file;
const char* method_name;
int use_checkpointing = 0;

// Parameter ranges and initial values
double param_lower_bound[NUM_PARAMS];
double param_upper_bound[NUM_PARAMS];
double param_current_value[NUM_PARAMS];

// Function to read calibration data from a CSV file
void read_calibration_data(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Could not open calibration data file");
        exit(EXIT_FAILURE);
    }

    char line[4096];
    int day = 0;

    // Skip header line
    fgets(line, sizeof(line), file);

    // Read data lines
    while (fgets(line, sizeof(line), file) && day < NUM_DAYS) {
        char *token = strtok(line, ",");
        // Skip date
        token = strtok(NULL, ",");

        // Read tidal range
        observed_data[day][0] = atof(token);
        token = strtok(NULL, ",");

        // Read salinity gradient for each cell
        for (int cell = 1; cell <= NUM_CELLS; cell++) {
            observed_data[day][cell] = atof(token);
            token = strtok(NULL, ",");
        }
        day++;
    }

    fclose(file);
}

// Function to read parameters from `params.txt`
void read_parameters(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Could not open parameter file");
        exit(EXIT_FAILURE);
    }

    char line[256];
    while (fgets(line, sizeof(line), file)) {
        if (sscanf(line, "use_checkpointing = %d", &use_checkpointing) == 1) continue;

        // Read parameter ranges and current values
        if (sscanf(line, "LC1 = %lf", &param_current_value[0]) == 1) continue;
        if (sscanf(line, "LC1_min = %lf", &param_lower_bound[0]) == 1) continue;
        if (sscanf(line, "LC1_max = %lf", &param_upper_bound[0]) == 1) continue;
        if (sscanf(line, "Chezy1 = %lf", &param_current_value[1]) == 1) continue;
        if (sscanf(line, "Chezy1_min = %lf", &param_lower_bound[1]) == 1) continue;
        if (sscanf(line, "Chezy1_max = %lf", &param_upper_bound[1]) == 1) continue;
        if (sscanf(line, "Rs1 = %lf", &param_current_value[2]) == 1) continue;
        if (sscanf(line, "Rs1_min = %lf", &param_lower_bound[2]) == 1) continue;
        if (sscanf(line, "Rs1_max = %lf", &param_upper_bound[2]) == 1) continue;

        if (sscanf(line, "LC2 = %lf", &param_current_value[3]) == 1) continue;
        if (sscanf(line, "LC2_min = %lf", &param_lower_bound[3]) == 1) continue;
        if (sscanf(line, "LC2_max = %lf", &param_upper_bound[3]) == 1) continue;
        if (sscanf(line, "Chezy2 = %lf", &param_current_value[4]) == 1) continue;
        if (sscanf(line, "Chezy2_min = %lf", &param_lower_bound[4]) == 1) continue;
        if (sscanf(line, "Chezy2_max = %lf", &param_upper_bound[4]) == 1) continue;
        if (sscanf(line, "Rs2 = %lf", &param_current_value[5]) == 1) continue;
        if (sscanf(line, "Rs2_min = %lf", &param_lower_bound[5]) == 1) continue;
        if (sscanf(line, "Rs2_max = %lf", &param_upper_bound[5]) == 1) continue;
    }
    fclose(file);
}

// Hypothetical objective function mimicking calibration (minimizing squared errors)
double objective_function(unsigned n, const double *x, double *grad, void *data) {
    iteration_count++;

    // Ensure values stay within reasonable bounds to avoid NaN
    for (unsigned i = 0; i < n; i++) {
        if (x[i] < param_lower_bound[i] || x[i] > param_upper_bound[i]) return HUGE_VAL;  // Prevent values outside bounds
    }

    double error = 0.0;

    // Simulate the model and compare with observed data
    for (int day = 0; day < NUM_DAYS; day++) {
        double simulated_tidal_range = 2.0 * x[0] + 0.5 * x[1] * sin((2 * M_PI * day) / NUM_DAYS);
        error += pow((simulated_tidal_range - observed_data[day][0]) / observed_data[day][0], 2);

        for (int cell = 1; cell <= NUM_CELLS; cell++) {
            double simulated_salinity = 30.0 + 5.0 * x[2] * sin((2 * M_PI * day) / NUM_DAYS) + 0.1 * x[3] * cos((2 * M_PI * cell) / NUM_CELLS);
            error += pow((simulated_salinity - observed_data[day][cell]) / observed_data[day][cell], 2);
        }
    }

    // Gradient calculation (if needed)
    if (grad) {
        // Placeholder for gradient calculation
        for (unsigned i = 0; i < n; i++) {
            grad[i] = 0.0; // Simplified; in practice, calculate the derivative of the error w.r.t. x[i]
        }
    }

    // Write parameter values to the output file
    fprintf(output_file, "%s,%d", method_name, iteration_count);
    for (unsigned i = 0; i < n; i++) {
        fprintf(output_file, ",%.10f", x[i]);
    }
    fprintf(output_file, ",%.10f\n", error);

    // Print progress every 100 iterations
    if (iteration_count % 100 == 0) {
        double progress = (double)iteration_count / max_iterations * 100;
        printf("🔄 Iteration %d | Progress: %.1f%% | Current error = %.8f\n", 
                iteration_count, progress, error);
    }

    // Stop early if max iterations reached
    if (iteration_count >= max_iterations) {
        printf("⚠️ Reached max iterations! Stopping optimization.\n");
        return HUGE_VAL;
    }

    return error;
}

void test_optimizer(nlopt_algorithm algorithm, const char *name) {
    iteration_count = 0;
    method_name = name;

    // Open the output file in append mode
    output_file = fopen("nlopt_log.csv", "a");
    if (!output_file) {
        perror("Could not open output file");
        exit(EXIT_FAILURE);
    }

    // Define the optimization problem
    nlopt_opt opt = nlopt_create(algorithm, NUM_PARAMS);
    if (!opt) {
        printf("❌ Failed to create NLopt optimizer for %s!\n", name);
        fclose(output_file);
        return;
    }

    // Define upper and lower bounds for parameters
    nlopt_set_lower_bounds(opt, param_lower_bound);
    nlopt_set_upper_bounds(opt, param_upper_bound);

    // Set the objective function
    if (nlopt_set_min_objective(opt, objective_function, NULL) < 0) {
        printf("❌ Failed to set objective function for %s!\n", name);
        nlopt_destroy(opt);
        fclose(output_file);
        return;
    }

    // Set optimization stopping criteria
    nlopt_set_xtol_rel(opt, 1e-6); // Higher precision for testing
    nlopt_set_maxeval(opt, max_iterations); // Prevent infinite loops

    // Initial guesses for parameters
    double x[NUM_PARAMS];
    for (int i = 0; i < NUM_PARAMS; i++) {
        x[i] = param_current_value[i];
    }
    double minf; // Minimum objective function value

    // Record start time
    clock_t start_time = clock();

    // Optimize
    printf("🔧 Starting optimization with %s...\n", name);
    nlopt_result result = nlopt_optimize(opt, x, &minf);

    // Record end time
    clock_t end_time = clock();
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    // Check the result
    if (result < 0) {
        printf("❌ NLopt optimization with %s failed! Error Code: %d\n", name, result);
    } else {
        printf("✅ NLopt found best parameters with %s and error = %.10f\n", name, minf);
        printf("   Optimal LC1 = %.2f | Chezy1 = %.2f | Rs1 = %.2f\n", x[0], x[1], x[2]);
        printf("   Optimal LC2 = %.2f | Chezy2 = %.2f | Rs2 = %.2f\n", x[3], x[4], x[5]);
        printf("   Time taken: %.2f seconds\n", elapsed_time);
    }

    // Write final results to the output file
    fprintf(output_file, "%s,Final", method_name);
    for (unsigned i = 0; i < NUM_PARAMS; i++) {
        fprintf(output_file, ",%.10f", x[i]);
    }
    fprintf(output_file, ",%.10f\n", minf);
    fprintf(output_file, "%s,Time,%.2f\n", method_name, elapsed_time);

    // Cleanup
    nlopt_destroy(opt);
    fclose(output_file);
}

int main() {
    printf("🔧 Starting NLopt Test...\n");

    // Open the output file and write the header
    output_file = fopen("nlopt_log.csv", "w");
    if (!output_file) {
        perror("Could not open output file");
        return 1;
    }
    fprintf(output_file, "Method,Iteration");
    for (int i = 0; i < NUM_PARAMS; i++) {
        fprintf(output_file, ",Param%d", i + 1);
    }
    fprintf(output_file, ",Error\n");
    fclose(output_file);

    // Read calibration data
    read_calibration_data("calibration_data.csv");
    printf("✅ Calibration data loaded.\n");

    // Read parameters
    read_parameters("INPUT/params.txt");
    printf("✅ Parameters loaded.\n");

    // Test different optimization algorithms
    test_optimizer(NLOPT_LN_COBYLA, "COBYLA");
    test_optimizer(NLOPT_LD_MMA, "MMA");
    test_optimizer(NLOPT_LN_SBPLX, "SBPLX");
    test_optimizer(NLOPT_LN_BOBYQA, "BOBYQA");

    printf("🔧 NLopt Test completed successfully.\n");
    return 0;
}