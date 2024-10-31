from Utils import Variable, Function
from loguru import logger
from pprint import pprint

class Var:
    r = Variable("r")
    t = Variable("t")
    x = Variable("x")
    y = Variable("y")
    z = Variable("z")
    w = Variable("w")
    s = Variable("s")
    H = Variable("H")
    xt = Variable("xt")
    yt = Variable("yt")
    tau1 = Variable("tau1")
    tau2 = Variable("tau2")
    alpha = Variable("alpha")
    beta = Variable("beta")
    gamma = Variable("gamma")

class Func:
    dxdt = Function("dxdt", "-1*r*x*y - s*x*z", [Var.r, Var.s, Var.x, Var.y, Var.z])
    dydt = Function("dydt", "-1*r*x*y + alpha*r*x*y*H", [Var.r, Var.alpha, Var.H, Var.x, Var.y])
    dwdt = Function("dwdt", "r*x*y", [Var.r, Var.x, Var.y])
    dzdt = Function("dzdt", "-1*x*s*z - gamma * z + beta * r * x * y * H", [Var.s, Var.gamma, Var.beta, Var.H, Var.x, Var.y, Var.z, Var.r])
    f1 = Function("f1", "x*y+w", [Var.x, Var.y, Var.w])
    f2 = Function("f2", "1.0 * 10e-12 + y + w", [Var.y, Var.w])
    
    dxdt_st1 = Function("dxdt_st1", "-r*x*y", [Var.r, Var.x, Var.y])
    dydt_st1 = Function("dydt_st1", "-r*x*y", [Var.r, Var.x, Var.y])
    dwdt_st1 = Function("dwdt_st1", "r*x*y", [Var.r, Var.x, Var.y])

    dxdt_st2 = dxdt_st1
    dydt_st2 = Function("dydt_st2", "-r*x*y + alpha*r*xt*yt*H", [Var.r, Var.alpha, Var.H, Var.xt, Var.yt])
    dwdt_st2 = dwdt_st1

class ModelInit: 
    x_zero = 0.75 * 10**(-4)
    y_zero = 10**(-15)
    w_zero = 0
    z_zero = 0
    alpha = 1.8
    beta = 20
    gamma = 0.002
    r = 5 * 10**4
    s = 1.0 * 10**5
    m1 = 0.3 * 10**(-13)
    m2 = 1.1 * 10**(-7)
    step = 0.1


class Model:
    def __init__(self, dxdt, dydt, dwdt, dzdt, f1, f2, x_zero, y_zero, w_zero, z_zero, alpha, beta, gamma, r, s, m1, m2, step):
        self.x_zero = x_zero
        self.y_zero = y_zero
        self.w_zero = w_zero
        self.z_zero = z_zero
        self.thau1 = 0
        self.thau2 = 0
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.r = r
        self.s = s
        self.m1 = m1
        self.m2 = m2
        self.H1 = 0
        self.H2 = 0
        self.step = step
        
        self.t0 = None
        self.t1 = None
        
        self.dxdt = dxdt
        self.dydt = dydt
        self.dwdt = dwdt
        self.dzdt = dzdt
        
        self.f1 = f1
        self.f2 = f2
        
        self.compute = {}
    
        
    
    
    @logger.catch
    def compute_first_state(self, mth=None):
        
        dxdt = Func.dxdt_st1
        dydt = Func.dydt_st1
        dwdt = Func.dwdt_st1
        step = self.step
        
        def euiler_step(self):
            x, y, w = self.x_zero, self.y_zero, self.w_zero
            self.x_zero += dxdt(self.r, x, y) * step
            self.y_zero += dydt(self.r, x, y) * step
            self.w_zero += dwdt(self.r, x, y) * step
        
        def runge_kutta_step(self,):
            x, y, w = self.x_zero, self.y_zero, self.w_zero
            k1x = dxdt(r=self.r, x=x, y=y) 
            k1y = dydt(r=self.r, x=x, y=y) 
            k1w = dwdt(r=self.r, x=x, y=y) 

            # Вычисление k2
            k2x = dxdt(r=self.r, x=x + 0.5 * step, y=y + 0.5 * step * k1y) 
            k2y = dydt(r=self.r, x=x + 0.5 * step, y=y + 0.5 * step * k1y) 
            k2w = dwdt(r=self.r, x=x + 0.5 * step, y=y + 0.5 * step * k1y) 

            # Вычисление k3
            k3x = dxdt(r=self.r, x=x + 0.5 * step, y=y + 0.5 * step * k2y) 
            k3y = dydt(r=self.r, x=x + 0.5 * step, y=y + 0.5 * step * k2y) 
            k3w = dwdt(r=self.r, x=x + 0.5 * step, y=y + 0.5 * step * k2y) 

            # Вычисление k4
            k4x = dxdt(r=self.r, x=x + 1.0 * step, y=y + 1.0 * step * k3y) 
            k4y = dydt(r=self.r, x=x + 1.0 * step, y=y + 1.0 * step * k3y) 
            k4w = dwdt(r=self.r, x=x + 1.0 * step, y=y + 1.0 * step * k3y) 
            
            # Обновление значений
            self.x_zero += (k1x + 2 * k2x + 2 * k3x + k4x) * (step / 6)
            self.y_zero += (k1y + 2 * k2y + 2 * k3y + k4y) * (step / 6)
            self.w_zero += (k1w + 2 * k2w + 2 * k3w + k4w) * (step / 6)
        
        integral = 0
        t0 = 0
        while integral < self.m1:
            self.compute[t0] = (self.x_zero, self.y_zero, self.w_zero)
            integral += self.f1(self.x_zero, self.y_zero, self.w_zero) * self.step
            if mth =='rk4':
                runge_kutta_step(self)
            else:
                euiler_step(self)
            t0 += step

        self.t0 = t0
        self.H1 = 1
        
    def compute_second_state(self):
        dxdt = Func.dxdt_st2
        dydt = Func.dydt_st2
        dwdt = Func.dwdt_st2
        step = self.step

        
   

    def __call__(self, mth = None):
        # Пример использования метода find_t0
        self.compute_first_state(mth)
        print("First phase passed.\nProliferation started at t0: ", self.t0)
        print('t0:', self.t0)
        self.compute_second_state()

# Создание модели
model = Model(Func.dxdt, Func.dydt, Func.dwdt, Func.dzdt, Func.f1, Func.f2, ModelInit.x_zero, ModelInit.y_zero, ModelInit.w_zero, ModelInit.z_zero, ModelInit.alpha, ModelInit.beta, ModelInit.gamma, ModelInit.r, ModelInit.s, ModelInit.m1, ModelInit.m2, ModelInit.step)

# Запуск модели
result = model()
pprint(model.compute)
