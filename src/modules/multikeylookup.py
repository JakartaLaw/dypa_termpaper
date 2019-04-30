class MultiKeyLookup():

    def __init__(self):
        self.memory = dict()

    def save(self, p, m, f, t, choice):
        key = (p, m, f, t)
        self.memory[key] = choice

    def lookup(self, p, m, f, t):
        key = (p, m, f, t)
        return self.memory[key]
