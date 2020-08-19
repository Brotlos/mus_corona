import simpy
import sim_epidemic as model

STEP = 1
END = 15


class EpidemicParameters(object):
    '''Parameters for epidemic calculation with respect to SIRD-model
        Epidemic-parameters:
            β (float): Infection rate
            γ (float): Recover rate
            δ (float): Death rate
            
            ωs (float): Quarantine rate for susceptible individuals
            ωi (float): Quarantine rate for infectious individuals
            ωe (float): Quarantine-End rate for susceptible individuals

            v (float): (natural) Birth rate
            μ (float): (natural) Death rate
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

    # v = 0.0 #Geburtenrate
    # μ = 0.0 #Allgemeine Sterberate
    # γ = 0.0 #Genesungsrate
    # κ = 0.0 #Rate potentiell infektiöser Kontakte
    # ωs = 0.0 #Quarantine rate for susceptible individuals
    # ωi = 0.0 #Quarantine rate for infectious individuals
    # ωe = 0.0 #Quarantine-end rate for susceptible individuals
    # q = 0.0 #Infektionswahrscheinlichkeit bei Kontakt
    # β = κ * q #None #Häufigkeit infektiöser Kontakte
    # # def get_β(self):
    # #     '''Calculate (β = κ*q) and return β if β == None otherwise returns β'''
    # #     return self.κ * self.q if not self.β else self.β
    # # λ = None #Infektionsrate #TODO aus parameters entfernen
    # # def get_λ(self, n):
    # #     '''Calculate (λ = β I/N) and return λ if λ == None otherwise returns β'''
    # #     return self.get_β() * self.i_class/n
    # δ = 0.0 #Sterberate (durch Infektion)

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
            return False
    def execute(self, population):
        population.params.set_parameters_from_event(self)
        if self.callback:
            self.callback(self, population)
        alive = self._keep_alive

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
    events = []

    #simulation results
    s_class_data = []
    i_class_data = []
    dn_class_data = []
    di_class_data = []
    r_class_data = []
    xs_class_data = []
    xi_class_data = []
    n_class_data = []

    def __init__(self, env, name, 
                 s_class_cap, i_class_cap = 0, d_class_cap = 0, r_class_cap = 0, xs_class_cap = 0, xi_class_cap = 0, 
                 epidemic_params = None, events = None):
        self.env = env

        self.name = name

        #setup params
        self.params = epidemic_params if epidemic_params else EpidemicParameters()#TODO defaults

        #event setup
        self.subscribe_event(events)

        #setup ressources
        self.s_class = s_class_cap#EpidemicClass(s_class_cap)#Number of susceptible individuals
        self.i_class = i_class_cap#EpidemicClass(i_class_cap)#Number of infectious individuals
        self.dn_class = d_class_cap#EpidemicClass(d_class_cap)#Number of naturally deceased individuals
        self.di_class = d_class_cap#EpidemicClass(d_class_cap)#Number of deceased infectious individuals
        self.r_class = r_class_cap#EpidemicClass(r_class_cap)#Number of recovered individuals
        self.xs_class = xs_class_cap#EpidemicClass(xs_class_cap)#Number of susceptible individuals in quarantine
        self.xi_class = xi_class_cap#EpidemicClass(xi_class_cap)#Number of infectious individuals in quarantine
        self.n_class = N = sum([self.s_class, self.i_class, self.r_class, self.xs_class, self.xi_class])
        
        self.s_class_data.append(self.s_class)
        self.i_class_data.append(self.i_class)
        self.dn_class_data.append(self.dn_class)
        self.di_class_data.append(self.i_class)
        self.r_class_data.append(self.r_class)
        self.xs_class_data.append(self.xs_class)
        self.xi_class_data.append(self.xi_class)
        self.n_class_data.append(self.n_class)

        # Start the run process everytime an instance is created.
        self.action = env.process(self.run())
    
    def subscribe_event(self, events):
        '''Append the given EpidemicEvent to the events-set of Population'''
        if not events is type(list):
            event = events 
            events = [event] #convert event in list with event
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
            
            #simulate model
            cS, cI, cR, cXs, cXi, cDn, cDi, cN = model.sirxd_update_population(S=self.s_class, I=self.i_class, R=self.r_class, 
                                        Xs=self.xs_class, Xi=self.xi_class, Dn=self.dn_class, Di=self.di_class, 
                                        β=self.params.β, γ=self.params.γ, δ=self.params.δ, κs=self.params.ωs, 
                                        κi=self.params.ωi, κe=self.params.ωe, v=self.params.v, μ=self.params.μ)
            
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


if __name__ == '__main__':
    #setup environment and populations
    env = simpy.Environment()
    params_ger = EpidemicParameters(β = 0.5, γ=0.03, δ= 0.01, ωs= 0.01, ωi= 0.1, ωe= 0.1,v= 0.0, μ = 0.005)
    pop_germany = Population(env, 'Deutschland', 80000000, epidemic_params= params_ger, events=EpidemicEvent('Event 1', 'env.now == 2', v=5))
    pop_germany.subscribe_event(EpidemicEvent('Event 2', 'env.now == 5', v=1))
    
    env.run(until=END)

    model.sirxd_plot(
        time=list(range(END+1)), 
        S=pop_germany.s_class_data,
        I=pop_germany.i_class_data, 
        R=pop_germany.r_class_data, 
        Xs=pop_germany.xs_class_data, 
        Xi=pop_germany.xi_class_data, 
        Dn=pop_germany.dn_class_data, 
        Di=pop_germany.di_class_data, 
        N=pop_germany.n_class_data)