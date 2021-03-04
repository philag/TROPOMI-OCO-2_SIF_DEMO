# TROPOMI/OCO-2 SIF Tutorial
The public github repository was created to support the [NASA ARSET training course](https://appliedsciences.nasa.gov/join-mission/training/english/arset-use-solar-induced-fluorescence-and-lidar-assess-vegetation) "Use of Solar Induced Fluorescence and LIDAR to Assess Vegetation Change and Vulnerability".

## Content
### Demo_presentation.jl
Is a [Pluto](https://github.com/fonsp/Pluto.jl) notebook to demonstrate how to read and select TROPOMI & OCO-2 SIF data for arbitrary spatial shapes, compute temporal averages, generate spatial composites (via oversampling), and evaluate measurement uncertainties.

### Case_Study_Illinois.jl
Is a [Pluto](https://github.com/fonsp/Pluto.jl) notebook to illustrate the impact of the 2019 flood in the US Midwest on the seasonality of photosynthesis over Illinois. 
## Instructions

1. Install [Julia](https://julialang.org/downloads/)

2. Clone this github directory: ```git clone https://github.com/philag/TROPOMI-OCO-2_SIF_DEMO```

3. Enter the project folder and start Julia with: ```julia --project=.```

4. Install all required packages and dependencies by: ```using Pkg; Pkg.instantiate()```

5. Start Pluto: ```using Pluto; Pluto.run(port=2345,launch_browser=false)```
   
   Now there should be link which can be opened with your favorite browser (Chrome works best)
   
   Note: If you are working remotely you need to open a ssh tunnel first: ```ssh username@host -N -L 2345:localhost:2345```