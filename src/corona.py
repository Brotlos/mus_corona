import simpy

STEP = 1
END = 15


class EpidemicClass(object):
    _capacity = int(0)
    def __init__(self, capacity = 0):
        self._capacity = capacity
    
    def put(self, put_amount):
        self._capacity += put_amount
    
    def get(self, get_amount):
        self._capacity -= get_amount
        #TODO Error
    
    def read(self):
        return self._capacity

class EpidemicParameters(object):
    '''Parameters for epidemic calculation with respect to SIRD-model'''
    def __init__(self, v=None, μ=None, γ=None, κ=None, q=None, δ=None, β=None, λ=None):
        self.v=v if v else self.v
        self.μ=μ if μ else self.μ
        self.γ=γ if γ else self.γ
        self.κ=κ if κ else self.κ
        self.q=q if q else self.q
        self.δ=δ if δ else self.δ
        self.β=β if β else self.β
        self.λ=λ if λ else self.λ

    v = 0.0 #Geburtenrate
    μ = 0.0 #Allgemeine Sterberate
    γ = 0.0 #Genesungsrate
    κ = 0.0 #Rate potentiell infektiöser Kontakte
    q = 0.0 #Infektionswahrscheinlichkeit bei Kontakt
    β=None #Häufigkeit infektiöser Kontakte
    def get_β(self):
        '''Calculate (β = κ*q) and return β if β == None otherwise returns β'''
        return self.κ * self.q if not self.β else self.β
    λ=None #Infektionsrate
    def get_λ(self, n):
        '''Calculate (λ = β I/N) and return λ if λ == None otherwise returns β'''
        return self.get_β() * self.i_class/n
    δ = 0.0 #Sterberate (durch Infektion)

    def set_parameters_from_event(self, event):
        '''Set the parameters based on the given Event'''
        self.set_parameters(v=event.v, μ=event.μ, γ=event.γ, κ=event.κ, q=event.q, δ=event.δ, β=event.β, λ=event.λ)

    def set_parameters(self, v=None, μ=None, γ=None, κ=None, q=None, δ=None, β=None, λ=None):
        '''Set the parameters based on the given values'''
        self.v=v if v else self.v
        self.μ=μ if μ else self.μ
        self.γ=γ if γ else self.γ
        self.κ=κ if κ else self.κ
        self.q=q if q else self.q
        self.δ=δ if δ else self.δ
        self.β=β if β else self.β
        self.λ=λ if λ else self.λ

class EpidemicEvent(object):
    """Provides parameters, name and condition. The condition can be based on 'env' or 'population'"""
    name=''
    condition=''
    alive = True
    def __init__(self, name, condition, keep_alive = False, v=None, μ=None, γ=None, κ=None, q=None, δ=None, β=None, λ=None, callback = None):
        self.v=v
        self.μ=μ
        self.γ=γ
        self.κ=κ
        self.q=q
        self.δ=δ
        self.β=β
        self.λ=λ
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
    '''Represents a population with class division according to the SIRD model'''
    events = []

    #simulation results
    s_class_data = []
    i_class_data = []
    d_class_data = []
    r_class_data = []

    def __init__(self, env, name, s_class_cap, i_class_cap = 0, d_class_cap = 0, r_class_cap = 0, epidemic_params = None, events = None):
        self.env = env

        self.name = name

        #setup params
        self.params = epidemic_params if epidemic_params else EpidemicParameters()#TODO defaults

        #event setup
        self.subscribe_event(events)

        #setup ressources
        self.s_class = EpidemicClass(s_class_cap)
        self.i_class = EpidemicClass(i_class_cap)
        self.d_class = EpidemicClass(d_class_cap)
        self.r_class = EpidemicClass(r_class_cap)

        # Start the run process everytime an instance is created.
        self.action = env.process(self.run())
    
    def get_N(self):
        '''Calculates the number of individuals in the population based on classes and returns the value'''
        return sum([self.s_class, self.i_class, self.d_class, self.r_class])
        
    def subscribe_event(self, events):
        '''Append the given EpidemicEvent to the events-set of Population'''
        if not events is type(list):
            event = events 
            events = [event]
        for event in events:
            if type(event) is EpidemicEvent:
                self.events.append(event)
            else:
                print('Event: {} konnte nicht abboniert werden'.format(event))

    def _execute_events(self):
        '''Check for events which should be triggered and execute'''
        #Alte einträge herausfiltern
        self.events = list(filter(lambda e: e.alive, self.events))
        current_events = list(filter(lambda e: e.check_condition(self), self.events))
        for event in current_events:
            event.execute(self)
            # self.params.set_parameters_from_event(event)

    def run(self):
        '''Population in process for SimPy'''
        while True:
            try:
                self._execute_events()
                
                print('An dieser Stelle kommt Kais Berrechnung. Zeit: {:3} v: {}'.format(self.env.now, self.params.v))

                if self.s_class.read() <= 0:
                    print("Die Population {} wurde ausgelöscht".format(self.name))
                    break
                yield env.timeout(STEP)
            except simpy.Interrupt:

                print('Population detected an Interrupt')



if __name__ == '__main__':
    #setup environment and populations
    env = simpy.Environment()
    pop_germany = Population(env, 'Deutschland', 80000000, events=EpidemicEvent('Event 1', 'env.now == 2', v=5))
    pop_germany.subscribe_event(EpidemicEvent('Event 2', 'env.now == 5', v=1))
    
    env.run(until=END)

     