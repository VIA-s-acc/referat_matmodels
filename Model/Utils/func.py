from .var import Variable

def remove_duplicates(lst):
    seen = set()
    result = []
    for item in lst:
        if item not in seen:
            result.append(item)
            seen.add(item)
    return result

class Function:
    def __init__[T: Variable](self, name, body, vars: T | list[T]) -> "Function":
        self.name = name
        self.body = body
        self.vars = remove_duplicates(vars) if isinstance(vars, list) else [vars]
    
    def __str__(self) -> str:
        return f"{self.name}({', '.join([str(var) for var in self.vars])}) = {self.body}"
    
    def __call__(self, *args, **kwargs) -> any:
        for var, arg in zip(self.vars, args):
            if var.name not in list(kwargs.keys()):
                kwargs[var.name] = arg
        try:
            return eval(self.body, kwargs)
        except Exception as e:
            print(f"Error in function {self}: {e} | check kwargs")
            return None
            
    def __add__(self, other: "Function"):
        if self.vars == other.vars:
            return Function(f'F', self.body +"+"+other.body, self.vars)
        else:
            raise ValueError("Cannot add functions with different variables")
    
    def __radd__(self, other: "Function"):
        return self.__add__(other)
    
    def __sub__(self, other: "Function"):
        if self.vars == other.vars:
            return Function(f'F', self.body +"-"+other.body, self.vars)
        else:
            raise ValueError("Cannot subtract functions with different variables")
        
    def __rsub__(self, other: "Function"):
        return self.__sub__(other)
    
    def __mul__(self, other: "Function"):
        if self.vars == other.vars:
            return Function(f"F", f"({self.body})" +"*"+f"({other.body})", self.vars)
        else:
            raise ValueError("Cannot multiply functions with different variables")
        
    def __rmul__(self, other: "Function"):
        return self.__mul__(other)
    
    def __truediv__(self, other: "Function"):
        if self.vars == other.vars:
            return Function(f'F', f"({self.body})" +"/"+f"({other.body})", self.vars)
        else:
            raise ValueError("Cannot divide functions with different variables")
        
    def __rtruediv__(self, other: "Function"):
        return self.__truediv__(other)
    
    def __floordiv__(self, other: "Function"):
        if self.vars == other.vars:
            return Function(f'F', f"({self.body})" +"//"+f"({other.body})", self.vars)
        else:
            raise ValueError("Cannot divide functions with different variables")
        
    def __rfloordiv__(self, other: "Function"):
        return self.__floordiv__(other)
    
    def __mod__(self, other: "Function"):
        if self.vars == other.vars:
            return Function(f'F', f"({self.body})" +"%"+f"({other.body})", self.vars)
        else:
            raise ValueError("Cannot divide functions with different variables")
        
    def __rmod__(self, other: "Function"):
        return self.__mod__(other)
    
    def __eq__(self, other: "Function") -> bool:
        return self.name == other.name and self.body == other.body and self.vars == other.vars
    
    def __repr__(self) -> str:
        return f"Function({self.name}, {self.body}, {self.vars})"
    
    def __hash__(self) -> int:
        return hash((self.name, self.body))
