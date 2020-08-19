import simpy


class Population (object):

    #region Parameter
    #TODO Parameter für Zeitabschnitt
    #endregion

    #TODO RESS

    def __init__(self, env, s_class_cap, i_class_cap = 0, r_class_cap = 0):
        self.env = env
        # Start the run process everytime an instance is created.
        self.action = env.process(self.run())
    
    def event_setup (self, events):
        #TODO übergeben events registrieren
        pass

    def run(self):
        while True:

            pass
            #Die Population existiert und infiziert sich gemäß der voreingestellten Parameter

            #TODO Die Ressourcen müssen nach den Parametern berrechnet werden.
            #Für jeden Zeitabschnitt

            #TODO Check for new events


def rate_changing_event(self, env, population):
    #TODO event auslösen für population
    pass



if __name__ == "__main__":
    env = simpy.Environment()
    pop_germany = Population(env, 80000000)

    pop_germany.run()