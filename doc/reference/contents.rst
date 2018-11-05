.. _UserGuide:

User Guide
==========

This is the user guide for the OmegaMaxEnt_TRIQS module. It should be sufficient ot obtain good results in many cases. For a more advanced use, which will sometimes be necessary, refer to the :math:`\Omega MaxEnt` `user guide`_.

Main function in OmegaMaxEnt_TRIQS
----------------------------------

Whether you have a scalar-, a matrix-valued or a block Green's function, the real frequency Green's function is obtained with the function **compute_GfReFreq()**, which takes an object of type Gf,GfImFreq, GfImTime, or a BlockGf containing objects of one of those types, and returns an object of type GfReFreq or a BlockGf containing GfReFreq objects. The function signature is::

    compute_GfReFreq(G, ERR=None, grid_params=[], name="$G^R(\omega)$", interactive_mode=True,
                     save_figures_data=True, save_G=True, comp_grid_params=[], non_uniform_grid=False,
                     inv_sym=False, mu=1, nu=1)

where the parameters are:

G:
    Gf, GfImFreq, GfImTime or BlockGf object.

    The input Matsubara Green function.

ERR:
    Optional numpy array.

    Standard deviation if G is scalar or has a single element.

    ERR must have the same shape as G.data and its indices must match the mesh.

    For a non-diagonal covariance, see the :math:`\Omega MaxEnt` `user guide`_.

grid_params:
    Optional list of the form :math:`[\omega_{min}, \Delta\omega, \omega_{max}]`.

    Defines the minimum frequency, the frequency step, and the maximum frequency, respectively, of the real frequency grid of the output Green function.

    If empty, the output grid is set by :math:`\Omega MaxEnt`.

    Note: this is only the output frequency grid. For the grid used in the calculation, see parameter comp_grid_params_.

name:
    Optional string.

    Name parameter of the returned GfReFreq object

interactive_mode:
    Optional boolean.

    Turns off the `interactive mode`_ of :math:`\Omega MaxEnt` if set to False.

save_figures_data:
    Optional boolean.

    Tells :math:`\Omega MaxEnt` not to save the figures files if set to False.

    Set to False if G is a matrix or a BlockGf.

save_G:
    Optional boolean.

    By default, the result is save in hdf5 format in file "G_Re_Freq.h5".

.. _comp_grid_params:

comp_grid_params:
    Optional list of the form [:math:`\Delta\omega`] or [:math:`\Delta\omega`, spectrum_width]` or [:math:`\Delta\omega`, spectrum_width, spectrum_center].

    Grid parameters used in the computation.

    :math:`\Delta\omega` is the frequency step used in the main spectral region, namely, the part of the grid where most of the spectral weight is located.

    spectrum_width is the width of the main spectral region (should be between 2 and 4 standard deviations typically).

    spectrum_center is the center of the main spectral region.

    :math:`\Delta\omega` and spectrum_width are ignored if not positive.

non_uniform_grid:
    Optional boolean.

    Tells :math:`\Omega MaxEnt` to use a non-uniform grid in the main spectral region for the computation.

    This accelerates the calculation if the spectrum has a peak at zero frequency that is very narrow compared to the total width of the spectrum.

inv_sym:
    Optional boolean.

    If G is a matrix or a BlockGf, set inv_sym to True if :math:`G_{ji}=G_{ij}`. This simplifies the calculation of the off diagonal elements.

mu, nu:
    Optional parameters involved in the calculation of off-diagonal elements of matrix-valued Green functions. See appendix C of the :math:`\Omega MaxEnt` `user guide`_ for more details.


:math:`\Omega MaxEnt` parameter files
--------------------------------------

:math:`\Omega MaxEnt` uses the file **OmegaMaxEnt_input_params.dat** to interact with the user. You can optionally create that file in advance with the function **create_params_file()**. This allows you to set some parameters that are not set by **compute_GfReFreq()**. Some parameters, which appear in the top section of the file, are exclusively set by **compute_GfReFreq()**. Otherwise you can modify the other parameters.

If **OmegaMaxEnt_input_params.dat** does not exist when **compute_GfReFreq()** is called, it will create it.

:math:`\Omega MaxEnt` also uses a file called **OmegaMaxEnt_other_params.dat**, also created by **create_params_file()**, which defines a certain number of parameters on which the computation depends. Do not modify this file unless you are a very advanced user.

All the parameters in **OmegaMaxEnt_input_params.dat** and **OmegaMaxEnt_other_params.dat** are described in the `user guide`_.

.. _`interactive mode`:


Interactive mode
----------------

In interactive mode, :math:`\Omega MaxEnt` displays figures during the execution. If parameter "display preprocessing figures" in **OmegaMaxEnt_input_params.dat** is enabled, figures are displayed during the preprocessing stage. Otherwise, figures are displayed only at the end of the calculation, showing the resulting Green function, along with different quantities useful to evaluate the quality of the results. Taking a look at those figures is important, especially when processing a set of data for the first time. You can thus ensure that the calculation is indeed complete, that the algorithm worked well, and assess the reliability of the result. Details about how to interpret the different quantities displayed are given in the user guide.

In interactive mode, the program pauses at the end of the preprocessing stage if parameter "preprocess only" is enabled, and at the end of the calculation. During the pause, you can modify **OmegaMaxEnt_input_params.dat** and resume the execution of :math:`\Omega MaxEnt` once the file is saved. Otherwise, if you are satisfied with the result displayed, you can exit the execution of :math:`\Omega MaxEnt` by closing all the figures and entering any character other than 'y' in the terminal to resume the execution of your python code. If interactive_mode=False, :math:`\Omega MaxEnt` will not display any figure and exit at the end of the calculation.

.. _`imaginary time`:

Imaginary time data
-------------------

:math:`\Omega MaxEnt` works internally with the Matsubara frequency Green function. Therefore, if the data are provided in imaginary time, the program first has to compute the Fourier transform of the Green function before starting the calculation. This calculation is fast because it is only one fast Fourier transform. However, if you provide errors, the covariance matrix must also be Fourier transformed. If the number of time points is a few hundreds at most, this calculation is also over quickly, but if the number of points is more than a thousand, the calculation time becomes of the order of a few minutes and more, which might seem long if you are working in interactive mode. Also, if the standard deviation does not depend on tau, the covariance matrix is also proportional to the identity matrix in Matsubara frequency. Therefore, do not provide any error in that case and no covariance matrix will be Fourier transformed. The result will not be affected because it does not depend on the absolute value of the standard deviation. On the other hand, if the standard deviation depends on tau, the covariance matrix has to be Fourier transformed. In case, and if you need to redo the calculation with the same data, you can accelerate the preprocessing by using the Fourier transform of your Green function that is saved as a GfImFreq object called 'G' in file "G_im_freq.h5" and the Fourier transform of the covariance matrix saved in files "covar_ReRe.dat", "covar_ImIm.dat" and "covar_ReIm.dat" in directory "Fourier_transformed_data". Use the lines "re-re covariance file", "im-im covariance file" and "re-im covariance file" in section INPUT FILES PARAMETERS of the file **OmegaMaxEnt_input_params.dat** to provide the covariance to :math:`\Omega MaxEnt`.


For more details on how to use :math:`\Omega MaxEnt`, see the user guide.


Example: Suppose you have saved a TRIQS Matsubara Green's function as 'G' in file "G.h5", here is a script to obtain the corresponding real frequency Green's function::


    from pytriqs.archive import HDFArchive as HA
    import OmegaMaxEnt_TRIQS as OT

    #load the Green's function
    A=HA("G.h5",'r'):
    G=A['G']

    GR=OT.compute_GfReFreq(G)


.. _OME_main_page: https://www.physique.usherbrooke.ca/MaxEnt/index.php/Main_Page
.. _`user guide`: https://www.physique.usherbrooke.ca/MaxEnt/index.php/User_Guide