from Utils import Variable, Function
import matplotlib.pyplot as plt
from loguru import logger
from pprint import pprint
from argparse import ArgumentParser
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
    wt = Variable("wt")
    tau1 = Variable("tau1")
    tau2 = Variable("tau2")
    tau1dt = Variable("tau1dt")
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

    tau1dt = Function("tau1dt", "(x+y+w)/(xt,yt,wt)", [Var.x, Var.y, Var.w, Var.xt, Var.yt, Var.wt])
    
    tau2dt = Function('tau2dt', "(1.0 * 10e-12 + y + w)/(1.0 * 10e-12 + yt + wt)", [Var.y, Var.w, Var.yt, Var.wt])

    dxdt_st2 = dxdt_st1
    dydt_st2 = Function("dydt_st2", "-r*x*y + alpha*r*xt*yt*H", [Var.r, Var.x, Var.y,Var.alpha, Var.H, Var.xt, Var.yt])
    dwdt_st2 = dwdt_st1
    

class ModelInit:
    x_zero = 5
    y_zero = 1
    w_zero = 0
    z_zero = 0
    alpha = 10
    beta = 20
    gamma = 0.2
    r = 0.02
    s = 0.05
    m1 = 15
    m2 = 15
    step = 0.1
    max_time = 100
    max_time_afterf3 = 50
    eps = 3e-2


class Model:
    def __init__(self, x_zero, y_zero, w_zero, z_zero, alpha, beta, gamma, r, s, m1, m2, step):
        self.x_zero = x_zero
        self.y_zero = y_zero
        self.w_zero = w_zero
        self.z_zero = z_zero
        self.tau1 = 0
        self.tau2 = 0
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
        
        self.xt1 = 0
        self.yt1 = 0
        self.xt2 = 0
        self.yt2 = 0
    
       
        self.int1 = None
        self.int2 = None

        self.t0 = None
        self.t1 = None
        
        self.stop = False
        self.T = 0

       
        self.compute = {}
   
       
   
   
    @logger.catch
    def compute_first_state(self):
       
        dxdt = Func.dxdt_st1
        dydt = Func.dydt_st1
        dwdt = Func.dwdt_st1
        step = self.step
       
        def euiler_step(self):
            x, y, w = self.x_zero, self.y_zero, self.w_zero
            self.x_zero += dxdt(self.r, x, y) * step
            self.y_zero += dydt(self.r, x, y) * step
            self.w_zero += dwdt(self.r, x, y) * step
       
        integral = 0
        integral2 = 0
        t0 = 0
        while integral < self.m1:
            self.compute[t0] = (self.x_zero, self.y_zero, self.w_zero, self.z_zero, 0, 0)
            integral += Func.f1(self.x_zero, self.y_zero, self.w_zero) * self.step
            integral2 += Func.f2(self.y_zero, self.w_zero) * self.step
            euiler_step(self)
            t0 += step
        
        self.int1 = integral
        self.int2 = integral2

        self.t0 = t0
        self.H1 = 1
     
    @logger.catch  
    def compute_second_state(self, mth=None):
        dxdt = Func.dxdt_st2
        dydt = Func.dydt_st2
        dwdt = Func.dwdt_st2
        tau1dt = Func.tau1dt
        step = self.step
        t1 = self.t0
        
        def euler_step(self):
            x, y, w, z, tau1, key = self.x_zero, self.y_zero, self.w_zero, self.z_zero, self.tau1, None
            if self.tau1 == 0:
                key = 0
            else:
                key = min(self.compute.keys(), key=lambda x:abs(x-self.tau1))
                    

            self.tau1 += Func.f1(x, y, z) / Func.f1(self.compute[key][0], self.compute[key][1], self.compute[key][2]) * step
            self.xt1 = self.compute[key][0]
            self.yt1 = self.compute[key][1]
            self.x_zero += dxdt(self.r, x, y) * step
            self.y_zero += dydt(self.r, x, y, self.alpha, self.H1, self.xt1, self.yt1) * step              
            self.w_zero += dwdt(self.r, x, y) * step

        integral = self.int2
        while integral < self.m2:
            self.compute[t1] = (self.x_zero, self.y_zero, self.w_zero, self.z_zero, self.tau1, 0)
            integral += Func.f2(self.y_zero, self.w_zero) * self.step
            
            if mth =='rk4':
                pass
            else:
                euler_step(self)
            
            t1 += step
            if t1 > ModelInit.max_time:
                self.compute[t1] = (self.x_zero, self.y_zero, self.w_zero, self.z_zero, self.tau1, 0)
                self.stop = True
                print('Max_time exceeded')
                break
            
        self.t1 = t1
        self.int2 = integral
        self.H2 = 1
        
        
    @logger.catch
    def compute_third_state(self):
        dxdt = Func.dxdt
        dydt = Func.dydt
        dwdt = Func.dwdt
        dzdt = Func.dzdt
        tau1dt = Func.tau1dt
        tau2dt = Func.tau2dt
        step = self.step

        t2 = self.t1
        def euler_step(self):
            x, y, w, z, tau1, tau2 = self.x_zero, self.y_zero, self.w_zero, self.z_zero, self.tau1, self.tau2
            if tau2 == 0:
                key2 = 0
            else:
                key2 = min(self.compute.keys(), key=lambda x:abs(x-self.tau2))
                
            key1 = min(self.compute.keys(), key=lambda x:abs(x-self.tau1))
            self.tau1 += Func.f1(self.x_zero, self.y_zero, self.z_zero) / Func.f1(self.compute[key1][0], self.compute[key1][1], self.compute[key1][2]) * step
            self.tau2 += Func.f2(self.y_zero, self.w_zero) / Func.f2(self.compute[key2][1], self.compute[key2][2]) * step
            self.xt1 = self.compute[key1][0]
            self.yt1 = self.compute[key1][1]
            self.wt1 = self.compute[key1][2]
            self.xt2 = self.compute[key2][0]
            self.yt2 = self.compute[key2][1]
            self.wt2 = self.compute[key2][2]
            
            self.x_zero += dxdt(self.r, self.s, x, y, z) * step
            self.y_zero += dydt(self.r, self.alpha, self.H2, x, y) * step
            self.w_zero += dwdt(self.r, x, y) * step
            self.z_zero += dzdt(self.s, self.gamma, self.beta, self.H2, x, y, z, self.r) * step
            
        while t2 < ModelInit.max_time_afterf3 + self.t1:
            self.compute[t2] = (self.x_zero, self.y_zero, self.w_zero, self.z_zero, self.tau1, self.tau2)
            euler_step(self)
            t2 += step
            if abs(self.x_zero - self.z_zero) < ModelInit.eps:
                self.T = t2
                
           
            
                
        self.compute[t2] = (self.x_zero, self.y_zero, self.w_zero, self.z_zero, self.tau1, self.tau2)
        

    def __call__(self):
        # Пример использования метода find_t0
        self.compute_first_state()
        print("First phase passed.\nProliferation started at t0: ", self.t0)
        self.compute_second_state()
        if self.stop == True:
            print("Second phase dont started.\nMax_time exceeded, t1: ", None)
            return
        print("Second phase passed.\nStarting plasmatic cells division at t1: ", self.t1)
        self.compute_third_state()
        print("Third phase passed.\nEnd of the model")
        print("T: ", self.T)

# Создание модели
model = Model(ModelInit.x_zero, ModelInit.y_zero, ModelInit.w_zero, ModelInit.z_zero, ModelInit.alpha, ModelInit.beta, ModelInit.gamma, ModelInit.r, ModelInit.s, ModelInit.m1, ModelInit.m2, ModelInit.step)
parser = ArgumentParser(prog='Waltman Model')
parser.add_argument('-p', nargs = '*', default = None, help = 'functions to plot', required=False)
args = parser.parse_args()


class Args:
    x,y,w,z,t1,t2 = False,False,False,False,False,False
    if args.p:
        x = 'x' in args.p
        y = 'y' in args.p
        w = 'w' in args.p
        z = 'z' in args.p
        t1 = 't1' in args.p
        t2 = 't2' in args.p

# Запуск модели
result = model()
# pprint(model.compute)
t = list(model.compute.keys())
x = [values[0] for values in model.compute.values()]
y = [values[1] for values in model.compute.values()]
w = [values[2] for values in model.compute.values()]
z = [values[3] for values in model.compute.values()]
tau1 = [values[4] for values in model.compute.values()]
tau2 = [values[5] for values in model.compute.values()]

plt.figure(figsize=(12, 6))
if Args.x:
	plt.plot(t, x, label='x(t) - concentration of free antigen molecules at time t ', color='blue')
if Args.y:
	plt.plot(t, y, label='y(t) - concentration of free receptor molecules at time t', color='orange')
if Args.w:
	plt.plot(t, w, label='w(t) - concentration of antigen bound to receptor molecules at time t ', color='green')
if Args.z:
	plt.plot(t, z, label='z(t) - concentration of free antibody molecules at time t ', color='yellow')
if Args.t1:
	plt.plot(t, tau1, label='tau1(t)', color='blue')
if Args.t2:
	plt.plot(t, tau2, label='tau2(t)', color='purple')




# Настройки графика

plt.title('Графики функций x(t), y(t), w(t), z(t)')
plt.xlabel('Время (t)')
plt.ylabel('Значения')
plt.legend()
plt.grid()
plt.show()

