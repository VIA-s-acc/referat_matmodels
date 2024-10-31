class Variable:
    def __init__(self, name):

        self.name = name
        
    def __str__(self):
        return self.name
    
    def __iter__(self):
        return iter([self])
    
    def __hash__(self):
        return hash(self.name)
    
    def __eq__(self, other):
        return self.name == other.name
    
    def __repr__(self):
        return self.name