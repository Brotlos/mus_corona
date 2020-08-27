from copy import deepcopy
import simpy

import sim_epidemic as model

STEP = 1 #Stepwidth; interpreted as days. Adjust the rates of EpidemicParameters to this value
END = 175*10 #Simulates this number of steps

class EpidemicParameters(object):
    '''Parameters for epidemic calculation with respect to SIRXD-model
        Epidemic-parameters:
            v (float): birth rate

            κ (float): Rate of potentially infectious contacts
            q (float): risk of infection at contact
            β (float):  = κ * q Infection rate; Frequency of infectious contacts
            γ (float): Recover rate
            δ (float): Death rate
            
            ωs (float): Quarantine rate for susceptible individuals
            ωi (float): Quarantine rate for infectious individuals
            ωe (float): Quarantine-End rate for susceptible individuals

            v (float): (natural) Birth rate
            μ (float): (natural) Death rate

            λ = β I/N: Force of infection
        Note:
            l/y corresponds to the duration of the infection
    '''

    def __init__(self, v=0.0, μ=0.0, γ=0.0, κ=0.0, ωs=0.0, ωi=0.0, ωe=0.0, q=0.0, δ=0.0, β=None, λ=None):
        self.v = v 
        self.μ = μ 
        self.γ = γ 
        self.κ = κ 
        self.ωs = ωs 
        self.ωi = ωi 
        self.ωe = ωe 
        self.q = q 
        self.δ = δ 
        self.β = β
        self.λ = λ

    def set_parameters_from_event(self, event):
        '''Set the parameters based on the given Event'''
        self.set_parameters(v=event.v, μ=event.μ, γ=event.γ, κ=event.κ, ωs=event.ωs, ωi=event.ωi, ωe=event.ωe, q=event.q, δ=event.δ, β=event.β, λ=event.λ)

    def set_parameters(self, v=None, μ=None, γ=None, κ=None, ωs=None, ωi=None, ωe=None, q=None, δ=None, β=None, λ=None):
        '''Set the parameters based on the given values'''
        self.v=v if v else self.v
        self.μ=μ if μ else self.μ
        self.γ=γ if γ else self.γ
        self.κ=κ if κ else self.κ
        self.ωs=ωs if ωs else self.ωs
        self.ωi=ωi if ωi else self.ωi 
        self.ωe=ωe if ωe else self.ωe
        self.q=q if q else self.q
        self.δ=δ if δ else self.δ
        self.β=β if β else self.β
        self.λ=λ if λ else self.λ

class EpidemicEvent(object):
    """Provides parameters, name and condition. The condition can be based on 'env' or 'population'"""
    name=''
    condition=''
    alive = True
    def __init__(self, name, condition, keep_alive = False, v=None, μ=None, γ=None, κ=None, ωs=None, ωi=None, ωe=None, q=None, δ=None, β=None, λ=None, callback = None):
        self.v = v
        self.μ = μ
        self.γ = γ
        self.κ = κ
        self.ωs = ωs 
        self.ωi = ωi 
        self.ωe = ωe 
        self.q = q
        self.δ = δ
        self.β = β
        self.λ = λ
        self.name = name
        self.condition = condition
        self._keep_alive = keep_alive
        self.callback = callback

    def check_condition(self, population):
        '''Evals the given condition'''
        try:
            env = population.env
            return eval(self.condition)
        except:
            print(f'Condition {self.condition} could not be evaluated')
            return False
    
    def execute(self, population):
        print(f'Event: {self.name} is executing at time: {population.env.now}')
        population.params.set_parameters_from_event(self)
        if self.callback:
            self.callback(self, population)
        self.alive = self._keep_alive

class Population (object):
    ''' Represents a population with class division according to the SIRXD model
        SIRXD-Classes:
            S/s_class (float): Number of susceptible individuals
            I/i_class (float): Number of infectious individuals
            R/r_class (float): Number of recovered individuals

            Xs/xs_class (float): Number of susceptible individuals in quarantine
            Xi/xi_class (float): Number of infectious individuals in quarantine
            Dn/dn_class (float): Number of naturally deceased individuals
            Di/di_class (float): Number of deceased infectious individuals

            -> with simulation data stored in <classname>_data
    '''

    def __init__(self, env, name, n_class_cap = None,
                 s_class_cap = None, i_class_cap = 1, dn_class_cap = 0, di_class_cap = 0, r_class_cap = 0,
                 xs_class_cap = 0, xi_class_cap = 0, epidemic_params = None, events = None):
        #region check and prepare __init__ args
        #At time zero there has to be min one individual in class infected
        if i_class_cap < 1: print('i_class_cap is to small: {}'.format(i_class_cap))
        #Alle rates has to be greater than zero
        if (n_class_cap < 0 or s_class_cap and s_class_cap < 0 or
            dn_class_cap < 0 or di_class_cap < 0 or
            r_class_cap < 0 or xs_class_cap < 0 or
            xi_class_cap < 0): print("Given args can't be zero or smaller")
        if n_class_cap:
            #if n_class_cap is set, this value lead with the other classes to s_class value
            s_class_cap = n_class_cap - sum([i_class_cap, r_class_cap,xs_class_cap, xi_class_cap])
        #endregion

        #set SimPy environment
        self.env = env

        #set meta
        self.name = name

        #setup params
        self.params = deepcopy(epidemic_params) if epidemic_params else EpidemicParameters()

        #event setup
        self.events = []
        self.subscribe_event(events)

        #setup ressources
        self.s_class = s_class_cap      #Number of susceptible individuals
        self.i_class = i_class_cap      #Number of infectious individuals
        self.dn_class = dn_class_cap    #Number of naturally deceased individuals
        self.di_class = di_class_cap    #Number of deceased infectious individuals
        self.r_class = r_class_cap      #Number of recovered individuals
        self.xs_class = xs_class_cap    #Number of susceptible individuals in quarantine
        self.xi_class = xi_class_cap    #Number of infectious individuals in quarantine
        self.n_class = N = sum([self.s_class, self.i_class, self.r_class, self.xs_class, self.xi_class])
        
        #simulation results
        self.s_class_data = [self.s_class]
        self.i_class_data = [self.i_class]
        self.dn_class_data = [self.dn_class]
        self.di_class_data = [self.di_class]
        self.r_class_data = [self.r_class]
        self.xs_class_data = [self.xs_class]
        self.xi_class_data = [self.xi_class]
        self.n_class_data = [self.n_class]

        # Start the run process everytime an instance is created.
        self.action = env.process(self.run())
    
    def subscribe_event(self, events):
        '''Append the given EpidemicEvent to the events-set of Population'''
        if not events is type(list) and type(events) == EpidemicEvent:
            event = events 
            events = [event] #convert event in list with event
        if not events:
            events = []
        for event in events:
            if type(event) is EpidemicEvent:
                self.events.append(event)
            else:
                print('Event: {} could not be subscribed'.format(event))

    def _execute_events(self):
        '''Check for events which should be triggered and execute'''
        self.events = list(filter(lambda e: e.alive, self.events))#filter entries which aren't alive

        current_events = list(filter(lambda e: e.check_condition(self), self.events))
        for event in current_events:
            event.execute(self)

    def run(self):
        '''Population in process for SimPy'''
        while True:
            self._execute_events()

            #auto adabt birth rate to deaths in infected classes: infected and infected in quarantine 
            # self.params.v = self.params.μ + (self.params.δ*(self.i_class+self.xi_class)*STEP)/self.n_class 
            
            #simulate model
            cS, cI, cR, cXs, cXi, cDn, cDi, cN = model.sirxd_update_population(S=self.s_class, I=self.i_class, R=self.r_class, 
                                        Xs=self.xs_class, Xi=self.xi_class, Dn=self.dn_class, Di=self.di_class, 
                                        β=self.params.β, γ=self.params.γ, δ=self.params.δ, κs=self.params.ωs, 
                                        κi=self.params.ωi, κe=self.params.ωe, v=self.params.v, μ=self.params.μ, dt=STEP)
            
            #save values of current step
            self.s_class = cS; self.s_class_data.append(cS)
            self.i_class = cI; self.i_class_data.append(cI)
            self.dn_class = cDn; self.dn_class_data.append(cDn)
            self.di_class = cDi; self.di_class_data.append(cDi)
            self.r_class = cR; self.r_class_data.append(cR)
            self.xs_class = cXs; self.xs_class_data.append(cXs)
            self.xi_class = cXi; self.xi_class_data.append(cXi)
            self.n_class = cN; self.n_class_data.append(cN)
            
            #checkup and next step
            if self.s_class <= 0:
                print("Die Population {} wurde ausgelöscht".format(self.name))
                break
            yield env.timeout(STEP)

def absolute_to_growrate(xn1, xn0, Δ01):
    return (xn1/xn0)**(1 / Δ01) - 1

def rate_relative(a, percentage):
    '''returns the rate relativ to another rate
    example: a = 0.1; b should be 4% of the total period
                => a will interpretted as 96% and b will be 0.1 / 96 * 4'''
    return a / (100 - percentage) * percentage

def plot_population(populations):
    '''Plot all given populations, each with its own figure.'''
    i=2; last = len(populations)
    for population in populations:
        model.sirxd_plot(
            time=list(range(END+1)), 
            title=population.name,
            S=population.s_class_data,
            I=population.i_class_data, 
            R=population.r_class_data, 
            Xs=population.xs_class_data, 
            Xi=population.xi_class_data, 
            Dn=population.dn_class_data, 
            Di=population.di_class_data, 
            N=population.n_class_data, figure=i, last_figure=last == i-1)
        i+=1

def set_β_with_R0(r0, γ):
    return r0 * γ

if __name__ == '__main__':

    #setup environment and populations
    env = simpy.Environment()

    #setup epidemic parameters
    params_ger = EpidemicParameters(β = 0.5, γ=0.03, δ= 0.01, ωs= 0.01, ωi= 0.1, ωe= 0.1, v= 0.01, μ = 0.005)
    params_ger_bare = EpidemicParameters(β = 0.5, γ=0.03, δ= 0.01, ωs= 0.0, ωi= 0.1, ωe= 0.1, v= 0.0051, μ = 0.005)
    
    #https://www.destatis.de/DE/Themen/Gesellschaft-Umwelt/Bevoelkerung/_inhalt.html
    pop_ger = 83.2 * 10 ** 6
    v_2019 = absolute_to_growrate(pop_ger + 778100, pop_ger, 365)
    μ_2019 = absolute_to_growrate(pop_ger + 939500, pop_ger, 365)
    β_cor = absolute_to_growrate(100000, 1, 30)
    γ_10days=1/10
    ωe_14days=1/14
    δ_08 = rate_relative(γ_10days, 9277/234853 * 100)#death rate of infection related to recover rate

    params_ger_2019 = EpidemicParameters(β=β_cor, γ=γ_10days, δ=δ_08, ωs=0, ωi=0, ωe=0, v= v_2019, μ = μ_2019)    
    params_ger_II = EpidemicParameters(β=β_cor, γ=γ_10days, δ= δ_08, ωs= 0.0, ωi= 0.0, ωe= 0.0, v= v_2019, μ = μ_2019)
    

    #setup populations
    pop_germany = Population(env, 'Deutschland',\
        n_class_cap = pop_ger, i_class_cap=1, epidemic_params= params_ger_2019)

    pop_germanyII = Population(env, 'Deutschland (mit Gegenmaßnahmen)', \
        n_class_cap = pop_ger, i_class_cap=1, epidemic_params= params_ger_2019)
    
    #events

    def end_lockdown(ev, population):
        now = population.env.now
        population.subscribe_event(EpidemicEvent('Ende Lockdown & Sensibilisierung & Quarantäne',\
        f'env.now == {now + 40}', β=β_cor * 0.75, ωi=0.1, ωs=0.0001, ωe=ωe_14days))#

    pop_germanyII.subscribe_event(EpidemicEvent('Lockdown & Quarantäne', callback=end_lockdown,\
        condition='population.i_class >= 1000000', β=set_β_with_R0(1,pop_germanyII.params.γ), ωi=0.1, ωs=0.001, ωe=ωe_14days))#pop_germany.params.β * 0.7
    

    # Start simulation
    print('Simulation started')
    env.run(until=END)
    print('Simulation finish succesfull')
    
    #plot
    plot_population([pop_germany, pop_germanyII,])

    print('Finish succesfull')