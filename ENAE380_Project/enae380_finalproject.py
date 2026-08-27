#The Planetary Transit Program!
#By Giselle Doh 
#University of Maryland College Park
#ENAE380 section 0104
from astropy import units as u #handles units of measurement
from astropy import time #handles astronomical time
from poliastro.bodies import Mercury,Venus,Earth, Mars, Jupiter,Saturn,Uranus,Neptune,Sun 
from poliastro.ephem import Ephem #position and velocities of the planets at different times
from poliastro.maneuver import Maneuver # for solving Lambert's problem
from poliastro.twobody import Orbit #for creating the orbit
from poliastro.util import time_range #for creating the timespan

from astropy.coordinates import solar_system_ephemeris #this is changing the ephemerides used to be from jpl data and not the defualt
solar_system_ephemeris.set("jpl")


class PlanetaryTransit:
    """Program for determining mission trajectories by solving Lambert's problem. Given a start and end planet,
    launch or arrival date, and a boundary date to complete the timespan, the program will calculate Lambert's problem 
    for every date in the timespan, then choose the date that corresponds to the lowest change in velocity. Features also include 
    determining the duration of mission using the best date for launch/arrival, and plotting the initial,transfer, and final orbit."""

    def __init__(self,starting_planet,arrival_planet,date_launch,date_arrival,boundary):
        self.start_planet=starting_planet
        self.end_planet=arrival_planet
        self.launch_date=self.get_time(date_launch)#format the end time
        self.arrival_date=self.get_time(date_arrival)#format the launch time
        self.boundary=self.get_time(boundary)#format the boundary time
        self.best_date=None
        self.timespan= self.create_timespan()
        self.orb0=None
        self.orbf=None
        print(f"Calculating your journey from {self.start_planet.name} to {self.end_planet.name}..." )
        pass

    def get_time(self,date):
        """A function to format the time string received by the user into an Astropy time 
        object,which is the format used to do all the Astropy and Poliastro caculations."""
        if date==None: 
            return None #return None if the original value for the date was None
        else:
            format_time=time.Time(date,scale='tdb')
            return format_time #if there is an actual date, format it 


    def create_timespan(self):
        """A function to create the timespan of potential arrival or launch dates"""
        if self.arrival_date==None: #if arrival date is none, timespan is from launch date to boundary
            time_span=time_range(self.launch_date,end=self.boundary,scale='tdb') 
        else: #if launch date is none, time span is from boundary to arrival date
            time_span=time_range(self.boundary,end=self.arrival_date,scale='tdb')

        return time_span        

    

    #GET TRAJECTORY CODE
    def get_trajectory(self):
        """Function for solving Lambert's problem for
        Multiple dates in the time span and select the date corresponding to the lowest delta v"""

        min_dv=None #will store value of the minimum dv
        maneuver=None #will store the final manuver object
       
        if self.arrival_date==None:
            self.orb0 = Orbit.from_ephem(Sun, Ephem.from_body(self.start_planet,self.launch_date), self.launch_date)

            for i in self.timespan:
                self.orbf = Orbit.from_ephem(Sun, Ephem.from_body(self.end_planet,i),i)
                #checking wether the time of flight is less than zero(would cause an error )
                tof = (self.orbf.epoch - self.orb0.epoch).to(u.s).value
                if tof <= 0:
                        continue  # skip invalid Lambert problems
                man_lambert = Maneuver.lambert(self.orb0, self.orbf)
                total_dv = man_lambert.get_total_cost()

                if min_dv==None or total_dv<min_dv: #check if current total dv<current min div
                    min_dv=total_dv #set the min=to the current total
                    self.best_date=i#set baset date to the corresponding date
                    maneuver=man_lambert#save manaeuver object for this best date and min dv
            print("The best date for arrival is",self.best_date.iso)        
        else:
            self.orbf = Orbit.from_ephem(Sun, Ephem.from_body(self.end_planet,self.arrival_date), self.arrival_date) 
            for i in self.timespan:#using a minimum findinder to look thrpugh dV for possible arrival/launch dates
                self.orb0 = Orbit.from_ephem(Sun, Ephem.from_body(self.start_planet,i),i) 
                tof = (self.orbf.epoch - self.orb0.epoch).to(u.s).value
                if tof <= 0:
                    continue  # skip invalid Lambert problems                    
                man_lambert = Maneuver.lambert(self.orb0, self.orbf)
                total_dv = man_lambert.get_total_cost()

                if min_dv==None or total_dv<min_dv:
                    min_dv=total_dv
                    self.best_date=i 
                    maneuver=man_lambert            
            print("The best date for launch is",self.best_date.iso)
        return maneuver

    def get_plot(self,man_lambert):
        """Function for plotting the initial, transfer, and final orbits"""
        from poliastro.plotting import OrbitPlotter3D #importing the class for 3d plotting
        op = OrbitPlotter3D()
        #applying the maneuver to the initial orbit, and retruning the transfer and final orbits
        orb_t, orbf = self.orb0.apply_maneuver(man_lambert, intermediate=True)
        op.plot(self.orb0, label="Initial orbit")
        op.plot(orb_t, label="Transfer orbit")
        op.plot(orbf, label="Final orbit")
        op._figure.show()
        return None

    def get_dt(self):
        """Function for finding the total duration of the mission
        after solving for the best date of launch/arrival. """
        # Created a dictionary to list the hours in one solar day on each planet. Source https://spaceplace.nasa.gov/days/en/
        #converted from hours(listed on the website) to seconds and formatted to Astropy time seconds
        times={
            Mercury:(1408*3600)*u.second,
            Venus:(5832*3600)*u.second,
            Earth:86400*u.second,
            Mars:(25*3600)*u.second,
            Jupiter:(10*3600)*u.second,
            Saturn:(11*3600)*u.second,
            Uranus:(17*3600)*u.second,
            Neptune:(16*3600)*u.second

        }
        if self.arrival_date ==None:#if arrival date is none, we found the best date for arrival
            self.arrival_date=self.best_date
        else: #if launch date is none, we found the best date for launch
            self.launch_date=self.best_date
        
        if self.end_planet in times:
            period=times[self.end_planet] #finding number of seconds in a day on the destination planet
        else:
            print("ERROR: Planet not found")
            exit()

        #calculating duration.Converting raw time object to seconds before converting to days
        duration=(self.arrival_date-self.launch_date).to(u.second)/period 
        
        #extracting raw values for the days and rounding to 2 decimal places
        print("The mission will take ",round(duration.value,2),"days in time on",self.end_planet.name,"!") 
       
        period=times[Earth]
        duration=(self.arrival_date-self.launch_date).to(u.second)/period
        print("The mission will take ", round(duration.value,2), "Earth days!")

def stop(user_input):
    """A function for checking if the user wants to exit the program.
    If the user enters a capital or lowercase e than the program will end."""
    
    if user_input.lower()=="e":
        print("Goodbye!")
        exit()
    return user_input    

def main():
    #creating a dictionary with the planets that matches them to specfifc inputs from the user
    planets={
        "0":Mercury,
        "1":Venus,
        "2":Earth,
        "3":Mars,
        "4":Jupiter,
        "5":Saturn,
        "6":Uranus,
        "7":Neptune

    }
    #asking for a starting planet
    start=stop(input("Choose a starting planet.\n" 
    "Options:\n" 
    "0 - Mercury\n" 
    "1 - Venus\n" 
    "2 - Earth\n" 
    "3 - Mars\n" 
    "4 - Jupiter\n" 
    "5 - Saturn\n" 
    "6 - Uranus\n" 
    "7 - Neptune\n"
    "type 'e' to quit\n"))
    #checking for correct input format
    if start in planets:
        start_planet=planets[start]
    else:
        print("ERROR:Planet not found")
        exit()
    
    #asking for destination planet
    dest=stop(input("Choose a destination.\n" 
    "Options:\n" 
    "0 - Mercury\n" 
    "1 - Venus\n" 
    "2 - Earth\n" 
    "3 - Mars\n" 
    "4 - Jupiter\n" 
    "5 - Saturn\n" 
    "6 - Uranus\n" 
    "7 - Neptune\n"
    "type 'e' to quit\n"))

#checking for correct planet format
    if dest in planets:
        end_planet=planets[dest]
        if end_planet==start_planet:
            print("ERROR: You cannot use the same planet for the start and end point")
            exit()

    else:
        print("ERROR: Planet not found")
        exit()

    #asking whether the user wants to provide a launch or arrival date
    aol=stop(input("Would you like to enter a Launch date or Arrival Date?\n" 
    "Options:\n" 
    "0 - Launch\n" 
    "1 - Arrival\n"
    "type 'e' to quit\n"))

    if aol=="0":
        launch_date=stop(input("Please enter the lanch date in the format 'YYYY-MM-DD'\n" 
        "type 'e' to quit\n"))
        arrival_date=None
    else:
        arrival_date=stop(input("Please enter the arrival date in the format'YYYY-MM-DD'\n " 
        "type 'e' to quit\n"))
        launch_date=None

   #asking user for the boundary date needed to create the timespan 
    boundary=stop(input("Please enter a boundary date to complete the range of possible\n " 
    "launch/arrival dates in the format 'YYYY-MM-DD'\n "
    "type 'e' to quit\n"))

    #checking to ensure that the timespan will be valid.
    if arrival_date==None:
        if boundary<=launch_date:
            print("ERROR: The boundary date can't be less than or equal to the launch date")
            exit()
    else:
        if boundary>=arrival_date:
            print("ERROR: The boundary date can't be greater than or equal to the arrival date")
            exit()   

    user=PlanetaryTransit(start_planet,end_planet,launch_date,arrival_date,boundary)

#getting the maneuver object to use for plotting
    man_lambert=user.get_trajectory()
    user.get_dt() #getting the duration of the full mission
    user.get_plot(man_lambert)#plotting the orbits


main()
