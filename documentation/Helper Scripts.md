# File: path_setup.m

* Called to add the paths to the workspace. Also adds the repository to the MATLAB path.
* Called from each of the figure scripts by using `run(fullfile(mfilename('fullpath'), '../../path_setup.m'))`
     * Adapt this line as needed to call from other directories

# File: makefreqmask.m

* path: `in_path/operational/general_funcs/makefreqmask.m`
* Creates a logical mask for the required frequencies