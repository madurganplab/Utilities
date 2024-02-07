# import Pkg; Pkg.add("LsqFit")
using LsqFit
using Plots

# Generate sample data from two exponential curves
x = range(0, stop=10, length=100)
y = 2 * exp.(-0.5 * x) + 0.5 * exp.(-2 * x) + 0.1 * randn(length(x))

# Define the model function as a sum of two exponential curves
function model_func(x, p)
    return p[1] * exp.(-p[2] * x) + p[3] * exp.(-p[4] * x)
end

# Define initial parameter values for the model function
initial_params = [1.0, 0.1, 1.0, 0.1]

# Perform least squares fitting
fit_result = curve_fit(model_func, x, y, initial_params)

# Extract estimated parameters and their uncertainties
estimated_params = fit_result.param
estimated_params_uncertainties = stderror(fit_result)
# Print estimated parameters and their uncertainties at 5 sigma deviations
sigma = 5
for i in 1:length(estimated_params)
    lower_bound = estimated_params[i] - sigma * estimated_params_uncertainties[i]
    upper_bound = estimated_params[i] + sigma * estimated_params_uncertainties[i]
    println("Estimated Parameter $i: ", estimated_params[i], " ± ", estimated_params_uncertainties[i], " at $sigma sigma deviations")
    println("Confidence Interval: [$lower_bound, $upper_bound]")
end

# Plot sample data
pp=scatter(x, y, label="Sample Data")

# Plot fitted function
x_fit = range(0, stop=10, length=100)
pp=plot!(x_fit, model_func(x_fit, estimated_params), linewidth=2, label="Fitted Function")

# Add legend
# legend("test")

# Show plot
display(pp)