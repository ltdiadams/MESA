# Modules For Experiments In Stellar Astrophysics: The Effects of Magnetic Fields/Magnetic Pressure on Massive Star Evolution

# Goals/Objectives:

Modify the “Modules for Experiments in Stellar Astrophysics” (MESA) stellar evolution code to incorporate an added pressure component introduced by the addition of a magnetic field in a massive rotating star. With this the definition of Pressure (P) in MESA becomes:

    P = P(radius) + P(gas) + P(magnetic)
    
We would like to define an initial magnetic field strength (B0) and have the code integrate this factor into the development of our star. The magnetic pressure should be calculated by the following equation:

    P(magnetic) = [(Magnetic Field)2] / [2∗ (permeability of space)] 
    P(magnetic) = B2 / 2∗ μ 0
    
Where μ 0 is a constant (https://www.mantaro.com/resources/impedance-calculator.htm), and B is found by applying some function of radius to the initial B0. (B = B0∗ f(r))

∗ All in CGS units∗ 


# Preparation:

While I cover some of the introduction here, it is a good idea to read the MESA papers: http://mesa.sourceforge.net/capabilities.html

It’s also likely that others have investigated some of what you’re trying to do before, you can search the mailing lists for questions people have asked which may be similar to yours at this link: https://lists.mesastar.org/pipermail/mesa-users/

There are a number of other resources but these are particularly useful in getting started with MESA:

https://jschwab.github.io/mesa-2017/

http://www.stellar-astrophysics.org/about.html

http://cococubed.asu.edu/mesa_market/guides.html

To get a good understanding of some of the physics involved you may want to read this paper: https://ui.adsabs.harvard.edu/abs/2009ARA% 26A..47..333D/abstract

# Code Capabilities and Use:

As of now for simplicity, we say that B = B0, so within MESA we can just use the initial surface magnetic field as defined within the inlist as “x_ctrl(1)” . You can use inlists in MESA to store initialization data for your stellar model: mass, rotation, metallicity, and other properties or parameters. In fortran these files would be called “namelists” .

This B value can be passed into another MESA file called “run_star_extras” , this is where the user inputs their own instructions to include physics in MESA. This file is particularly powerful given its ability to modify input physics using add-ons known as “other” hooks. You can change these routines to account for the added pressure term for example.

To call this value, write the following fortran code or something similar:

    B = s% x_ctrl(1)
    
The ‘s’ is a pointer to your star, a “star pointer” . You can then use B0 within the subroutine
where it was defined.
Pressure is defined in a private file labeled “micro.f90” , while a more permanent addition of the pressure term would involve reforms to the module which includes this, for our purposes, making use of “run_star_extras” should do fine.
In the case of this project, what was attempted was addition of magnetic pressure to the overall pressure within the “other_mlt” hook. Done by including: 

     s% P(k) = s% Prad(k) + s% Pgas(k) + P_mag or similarly, 
     s% P(k) = s% P(k) + P_mag, where k represents “zones” within the star.
     
This pressure can also be used in of itself as a new term within the “data_for_extra_profile_columns” subroutine.
The variables you can use in “run_star_extras” are defined in “star_data.inc” with units. The directory “LOGS” within your working directory is very important because this is where
MESA output is stored when it is run.

Within LOGS you can notice that there are two varieties of files, the history.data file as well as the profile# .data files.

The “history.data” file contains information on the time-evolution of global model parameters.

The “profile# .data” files contain information on the internal structure of the star at certain times.

If you want more data to be outputted to the history.data file, this can be done with a few lines of code in the data_for_extra_history_columns subroutine. And with a small change to the how_many_extra_history_columns subroutine as well.
In the subroutine “how_many_extra_history_columns” the line how_many_extra_history_columns = # , has to be changed to the number of columns you wish to add, currently there is only
1 being added to MESA.

Then in the data_for_extra_history_columns subroutine, the following was added.

        B0 = s% x_ctrl(1)
        mu_not = (4∗ pi∗ 1.e-7)∗ ((1.e7)/(4∗ pi)) 
        P_mag = (B0∗ ∗ 2)/(2∗ mu_not) names(1) = ’P_mag’
        vals(1) = P_mag
        
∗ Variables should be initialized at the top of the routine∗

Similarly, if you want to add a term into the profile files; it requires the same change to the subroutine “how_many_extra_profile_columns” .
Then in “data_for_extra_profile_columns” something like the following could be included:

        if (ierr /= 0) return names(1) = "P"
        call star_ptr(id, s, ierr) if (ierr /= 0) return
        B0 = s% x_ctrl(1)
        mu_not = (4∗ pi∗ 1.e-7)∗ ((1.e7)/(4∗ pi)) P_mag = (B0∗ ∗ 2)/(2∗ mu_not)
        do k=1, s% nz                                                                 ! nz means number of zones
           P_with = s% P(k) vals(k, 1) = P_with
        enddo
        
Now, “P” will be included in the profile files.


MESA Reader is a tool which is particularly useful in making the plots that we’ve come up with thus far in the project. The “readmesa.py” script was used to plot data from the profile files in terms of the radius of the star in solar radii. (ie. log_opacity, log_conv_vel, v_div_csound, etc. all vs. stellar radius). This is better for our purposes than the pgstar functionality built into MESA.

The initialization of the star for these plots was a mass of 15 Msun (when M > 7-8 Msun star is “massive” ). And the value for B0 was 1000 Gauss to be used in calculating magnetic pressure.

In order to use MESA Reader, you can first open the “profile_columns.list” file and uncom- ment all the variables which you want MESA to include in the profile.data files. They can also be found at this link: http://www.ster.kuleuven.be/∼ pieterd/mesa/profile_columns.html.

These variables can be used in the script for plotting very simply! Just be sure to make use of python’s reload function in iPython when plotting multiple lines.


# Important Commands:

    ipython –pylab
    from importlib import reload
    import readmesa

* make sure you’re in the directory of readmesa.py

* cd to working directory

      readmesa.make_HR("LOGS/history.data", 1)
    
      readmesa.read_profile(“LOGS” , “profile# .data” , “label-string” )


Note: the label functionality is an alteration to the initial readmesa script.
More information on MESA Reader can be found here: https://billwolf.space/py_mesa_reader/MesaData.html

# Further exploration:

As of now, there are still problems which need to be resolved. When magnetic pressure is turned on, the model converges but does not run to completion, only to about profile54. And while MESA will run the standard model (no added pressure) all the way through, radius expands when it should not. With this, the plots seem to show some noise which should not exist meaning there could be some issues with the incorporation of the magnetic pressure in MESA.

Once the problems have been resolved we would like to add an inclination to the dipole magnetic field (oblique dipole) as well as consider other stellar properties, like: thermody- namic consistencies, perhaps impact of rotation on stability of magnetic fields, how fossil fields actually evolve with a star, etc.
