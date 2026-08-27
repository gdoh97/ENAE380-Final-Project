ENAE380 Final Project: Planetary Transfer Program

Project Overview:

My final project is a Python program that will calculate the best spacecraft launch/arrival date based on the user inputs by solving Lambert’s problem for a range of boundary dates. The program also tells the user the duration of the entire mission using this best date in Earth days and in days on the destination planet. Finally, the code creates a 3D plot of the initial orbit, transfer orbit, and final orbit and displays them in the default browser of the user.

User Inputs:
- Departure planet (int: 0 = Mercury, 1 = Venus, 2 = Earth, 3 = Mars, 4 = Jupiter, 5= Saturn, 6 = Uranus, 7 = Neptune)
- Arrival planet (int: 0 = Mercury, 1 = Venus, 2 = Earth, 3 = Mars, 4 = Jupiter, 5= Saturn, 6 = Uranus, 7 = Neptune)
- Earliest departure date (YYYY-MM-DD)
- Latest arrival date (YYYY-MM-DD)

Project Structure:
- enae380_finalproject.py → Project Program
- environment.yml → conda environment 

How to Install and Run:
1. Create the conda environment using the included environment.yml file:

```
conda env create -f environment.yml
conda activate enae380
```

2. Activate the conda environment and run the script:
```
python enae380_finalproject.py
```


Dependencies:
- Poliastro
- Astropy
- Matplotlib
- NumPy

