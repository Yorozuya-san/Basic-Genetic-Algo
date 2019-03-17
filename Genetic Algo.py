import numpy as np
from random import randint
from numpy.random import normal

class GeneticAlgo:

    def __init__(self, pop_size, number_bits, number_of_gen=10, prob_mutation=0.8):
        self.number_of_gen = number_of_gen
        self.pop_size = pop_size
        self.number_bits = number_bits
        self.prob_mutation = prob_mutation

        self.init_pop()
        self.execution()

    def conver_bin(self,a):
        x = []
        while a>0:
            x.append(a%2)
            a //= 2
        x.reverse()
        return x

    def init_pop(self):
        #taking samples over the distributtion
        #such that no. of 0 almost= no. of 1
        
        tot = int(2**self.number_bits)
        u = tot // self.pop_size
        l = 1
        self.pop = []
        for i in range(self.pop_size):
            p =[]
            num = randint(l,u)
            p = self.conver_bin(num)
            p = ([0]*(self.number_bits-len(p))) + p
            l = u+1
            u = (i+2) * (tot//self.pop_size)
            self.pop.append(p)

    def deci(self,x):
        a = 0
        number_bits = len(x)
        for i in range(number_bits):
            a += x[i] * (2**(number_bits-i-1))
        return a

    def fitness(self,chromo):
      ######### Ur one line fitness function goes here ########
      decimal_value = self.deci(chromo)
      x = decimal_value
      
      val = None #function
      
      #example
      #val = x**2  
      ################################################

      return val

    def find_index_value(self,x,prob):
        for i in range(len(prob)):
            if x<prob[i]:
                return i
        return len(prob)-1

    def roul(self,fit):
        tot = sum(fit)
        prob = [round((fit[0]/tot)*1000)]
        
        for i in range(1,self.pop_size):
            x = fit[i]/tot
            x = round(x*1000)
            prob.append(prob[-1]+x)
            
        sel = []
        len_prob = len(prob)
       
        #using the *1000 and cumilative freq one
        for i in range(self.pop_size//2):
            x = self.find_index_value(randint(1,999),prob)
            sel.append(self.pop[x])
            y = self.find_index_value(randint(1,999),prob)
            while x==y:
                y = self.find_index_value(randint(1,999),prob)
            sel.append(self.pop[y])
        
        return sel

    def cross_over(self, selected):
        cross_over_sites = []
        for i in range(self.pop_size//2):
            x = randint(1,self.number_bits-2)
            while x in cross_over_sites:
                x = randint(1,self.number_bits-2)
            cross_over_sites.append(x)
            
        self.pop = []
        
        for i in range(self.pop_size//2):
            a = selected[i*2][:cross_over_sites[i]] + selected[i*2+1][cross_over_sites[i]:]
            b = selected[i*2+1][:cross_over_sites[i]] + selected[i*2][cross_over_sites[i]:]
            self.pop.append(a)
            self.pop.append(b)

    def mutation(self):
        tot = 0
        for i in self.pop:
            for j in range(len(i)):
                if normal(3,1) < self.prob_mutation:
                    i[j] = (i[j]+1)%2
                    tot += 1
        return tot

    def analysis(self):
        # example of analysis
        ''' 

        pop_deci = []
        for i in self.pop:
            pop_deci.append(self.deci(i))
            
        return sum(pop_deci)/self.pop_size,max(pop_deci) #returning average and maximum of population
        '''
        

    def execution(self):
        for i in range(self.number_of_gen):
            fit = []
            
            #============== eval fitness ======================
            for j in range(self.pop_size):
                fit.append(self.fitness(self.pop[j]))
            
            #============== selecting the chormosomes using roulette wheel method ======================
            selected = self.roul(fit)
            
            print('\n\Generation',i)
            print('\nSelected pop',selected)
            
            #============== crossover of chromosomes ======================
            self.cross_over(selected)
            
            print('\nNew pop', self.pop)
            
            #============== mutation ====================================
            tot_mu = self.mutation()
            print('\nAfter mutation', self.pop)
            print('\nTotal mutations', tot_mu)
            
            #============== analysis of new generation ======================
            avg, mx = self.analysis()
            
            print('\nMax : {}\nAvg : {}'.format(mx,avg))
        
            

ga = GeneticAlgo(number_of_gen=10,pop_size=4,number_bits=6)



