import simpy
import numpy as np
import matplotlib.pyplot as plt
from enum import Enum
from random import random, randint, seed, choices, gauss

# init random number generator 
# to get reproducible results
seed(42)

# Print debugging output 
debugging = False
def debug(*args, **kwargs):
    global debugging
    if debugging:
        print(*args, **kwargs)

def randbool(probability):
    return True if random() < probability else False

def randlocation(x_min=0, x_max=50, y_min=0, y_max=50):
    x, y = randint(x_min, x_max), randint(y_min, y_max)
    return (x, y)

def distance_squared(location1, location2):
    diff = [v1-v2 for v1, v2 in zip(location1, location2)]
    dist = sum(v**2 for v in diff)
    return dist

# enum
class SIR(Enum):
    susceptible = 0
    infectious = 1
    recovered = 2


# simulation parameters
mu_old_age = 24 * 7
sigma_old_age = 24 * 2
people = []
num_people = 100 * 5
initial_num_infectious = 3
infectious_distance_squared = 8**2
infection_probability = 0.1
mu_infection_duration = 24 * 3
sigma_infection_duration = 24 * 1
# vaccination_probability = 0.4
sim_time = 1000
groups_sample_time = 1
stats_sample_time = 24


T, S, I, R = [], [], [], []
β, λ, γ, R0, Reff = [0], [0], [0], [0], [0]
def update_groups(env):
    global people
    global T, S, I, R 
    global β, λ, γ, R0, Reff
    global groups_sample_time

    while True:
        debug("Updating SIR Groups")
        cS, cI, cR = 0, 0, 0
        # N = len(people)
        for person in people:
            if person.state == SIR.susceptible:
                cS += 1
            elif person.state == SIR.infectious:
                cI += 1
            elif person.state == SIR.recovered:
                cR += 1

        T.append(env.now)
        S.append(cS)
        I.append(cI)
        R.append(cR)

        if cI == 0:
            # terminate life processes
            for person in people:
                person.life_process.interrupt("No more infectious")
            break

        yield env.timeout(groups_sample_time)

def update_stats(env):
    global num_people
    global β, λ, γ, R0, Reff
    global stats_sample_time

    yield env.timeout(stats_sample_time)
    while True:

        try:

            debug("Updating Stats")
            cS, cI, cR = S[-1], I[-1], R[-1]
            dS, dI, dR = (cS - S[-(stats_sample_time+1)]), (cI - I[-(stats_sample_time+1)]), (cR - R[-(stats_sample_time+1)])

            cλ = -dS / cS if cS else 0
            cβ = cλ * num_people / cI if cI else 0
            cγ = dR / cI if cI else 0
            cR0 = cβ / cγ if cγ else 0
            cReff = cR0 * cS / num_people

            λ.append(cλ)
            β.append(cβ)
            γ.append(cγ)
            R0.append(cR0)
            Reff.append(cReff)


            if cI == 0:
                break

        except IndexError:
            pass

        yield env.timeout(stats_sample_time)

def num_infectious():
    global I
    return I[-1]

class Person(object):

    @classmethod
    def gen_name(cls):
        # global people_counter
        try: cls.people_counter
        except: cls.people_counter = 0

        cls.people_counter += 1
        return cls.people_counter

    @classmethod
    def gen_home(cls):
        return randlocation()

    def __init__(self, env, name=None, home=None, state=SIR.susceptible):
        global mu_old_age, sigma_old_age, mu_infection_duration, sigma_infection_duration

        self.env = env
        self.name = name if name else self.gen_name()
        self.home = home if home else self.gen_home()
        self.location = self.home
        self.is_outside = False
        
        self.birth_time = env.now
        self.lifespan = gauss(mu_old_age, sigma_old_age)

        self.state = state
        self.infection_time = 0
        self.recover_time = max(gauss(mu_infection_duration, sigma_infection_duration), 24)

        # self.in_quarantine = False
        # self.is_vaccinated = False

    def __call__(self):
        return self.life()

    def age(self):
        return (self.env.now - self.birth_time)

    def life(self):
        debug(f"P{self.name} is born.")

        try:
            while True:

                # Determine daily routine
                sleep_time = randint(4, 8)
                day_time = randint(4, 8)
                is_going_outside = bool(random() <= 0.5)

                # aging
                if self.age() >= self.lifespan:
                    debug(f"P{self.name} died.")
                    self.rebirth()

                # sleep
                yield self.env.process(self.sleep(sleep_time))

                # go outside or stay at home
                if is_going_outside:
                    # Go Outside
                    outside_location = randlocation()
                    yield self.env.process(self.go_outside(day_time, outside_location))
                else:
                    # Stay Home
                    yield self.env.process(self.stay_home(day_time))

        except simpy.Interrupt:
            return

    def sleep(self, sleep_time):
        debug(f"P{self.name} goes to sleep for {sleep_time}.")
        yield self.env.timeout(sleep_time)
        self.recover()
        debug(f"P{self.name} woke up at {self.env.now} after {sleep_time}.")

    def stay_home(self, home_time):
        debug(f"P{self.name} stays home for {home_time}.")
        yield self.env.timeout(home_time)
        debug(f"P{self.name} stayed home for {home_time} until {self.env.now}.")

    def go_outside(self, outside_time, outside_location):
        global people, infectious_distance_squared, infection_probability
        # go to location
        self.location = outside_location
        self.is_outside = True
        debug(f"P{self.name} goes outside to {outside_location}.")

        # infecting yourself
        if self.state == SIR.susceptible:
            # find out if infectious people are nearby
            local_infectious = [person for person in people \
                if person.is_outside and person.state == SIR.infectious and \
                distance_squared(self.location, person.location) <= infectious_distance_squared]

            # if infectious people are nearby, get infected
            if local_infectious and randbool(infection_probability):
                self.get_infected()
        # infecting others
        elif self.state == SIR.infectious:
            # find out if susceptible people are nearby
            local_susceptible = [person for person in people \
                if person.is_outside and person.state == SIR.susceptible and \
                distance_squared(self.location, person.location) <= infectious_distance_squared]

            for person in local_susceptible:
                # if infectious people are nearby, get infected
                if randbool(infection_probability):
                    person.get_infected()

        # stay outside for some time
        yield self.env.timeout(outside_time)

        # go home
        debug(f"P{self.name} goes home to {self.home}.")
        self.is_outside = False
        self.location = self.home
    
    def rebirth(self):
        self.__init__(self.env)
        debug(f"P{self.name} is born.")


    def get_infected(self): 
        self.state = SIR.infectious
        self.infection_time = self.env.now
        debug(f"P{self.name} got infected.")       

    def recover(self):
        if self.state == SIR.infectious:
            infectious_time = self.env.now - self.infection_time
            if infectious_time >= self.recover_time:
                self.state = SIR.recovered
                debug(f"P{self.name} recovered.")

def plot_results():
    global T, S, I, R 
    global β, λ, γ, R0, Reff
    print(f"Plotting {len(T)} data points.")

    plt.figure()
    plt.stackplot(T, R, I, S, labels=("Recovered, Infectious, Susceptible"), colors=("green", "orange", "blue"))
    plt.title("SIR Stackplot")
    plt.xlabel("Time [h]")
    plt.ylabel("People")
    plt.xlim(0, len(T))
    plt.legend((
        "Recovered", 
        "Infectious", 
        "Susceptible", 
        # "Susceptible in quarantine", 
        # "Infectious in quarantine", 
        # "Naturally deceased", 
        # "Deceased infectious", 
        # "Total population"
    ))

    plt.figure()
    plt.plot(T, S, T, I, T, R)
    plt.title("SIR Plot")
    plt.xlabel("Time [h]")
    plt.ylabel("People")
    plt.xlim(0, len(T))
    plt.legend((
        "Susceptible", 
        "Infectious", 
        "Recovered", 
        # "Susceptible in quarantine", 
        # "Infectious in quarantine", 
        # "Naturally deceased", 
        # "Deceased infectious", 
        # "Total population"
    ))

    stats_T = range(len(R0))
    plt.figure()
    plt.plot(stats_T, β, stats_T, λ, stats_T, γ)
    plt.title("SIR β, λ, γ")
    plt.xlabel("Time [d]")
    plt.ylabel("β, λ, γ")
    plt.xlim(0, len(stats_T))
    plt.legend((
        "β", 
        "λ", 
        "γ", 
    ))

    plt.figure()
    plt.plot(stats_T, R0, stats_T, Reff)
    plt.title("SIR R0, Reff")
    plt.xlabel("Time [d]")
    plt.ylabel("R0, Reff")
    plt.xlim(0, len(stats_T))
    plt.legend((
        "R0", 
        "Reff", 
    ))
    plt.show()


if __name__ == "__main__":

    # Create simulation environment
    env = simpy.Environment()

    # Processes that update SIR & stats
    groups_process = env.process(update_groups(env))
    stats_process = env.process(update_stats(env))

    # Create Person
    people = [Person(env) for i in range(num_people)]

    # Add Person to simulation environment
    for person in people:
        person.life_process = person.env.process(person.life())

    # initial infectious
    for person in choices(people, k=initial_num_infectious):
        person.get_infected()

    # rum Simulation
    print("Running simulation")
    print("This might take some time...")
    env.run(sim_time)

    # plot simulation results
    plot_results()