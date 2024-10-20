# SUN_SUT_Plotter
Tools implemented for the visualization of SUN and SUT distributions


The SUN_SUT_Plotter provides a efficient way to visualize SUN and SUT distribution. Before using the code, please make sure that the Psi and Gamma in the dp list should have capitals for the first letter. x_lim and y_lim can be adjusted in the input parameters. You can specify one or both, controlling the boundaries of the Y_1 and Y_2. If x_lim and y_lim have no input values, they will be set to the default boundaries. I have implemented two functions in the end A_to_H and H_to_A. They do the change of parameterizations by taking a list of parameters dp. In addition, make sure you run the dsut and rsut functions first.
