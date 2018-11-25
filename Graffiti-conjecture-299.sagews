︠11c4573a-d5a1-4053-bbcc-fdbd3a9a20ba︠
︡2e03b1f6-6731-48ba-bd65-fd2e1571e0be︡
︠e5d5daf5-7a3b-40e0-9552-1786b4687236s︠
# Conjecture: total_domination_number (G) <= annihilation(G) + powerset_S(G)

def degrees(G):
    stopnje = G.degree_sequence()
    return stopnje

def powerset_S(G):
    velikost_S = degrees(G).count(2)
    return velikost_S


def annihilation(G):
    stopnje = degrees(G)
    stopnje_obr = stopnje[::-1]
    vel = int(sum(stopnje_obr)/2)
    sest = 0
    if stopnje_obr == [0]: #če ima graf samo eno vozlišče
        return 1
    else:
        for i, k in enumerate(stopnje_obr):
            sest += k
            if sest > vel:
                return i


def total_domination_number(G):
    if len(G.vertices()) == 1:
        return 0
    elif len(G.vertices()) == 2:
        return 1
    p = MixedIntegerLinearProgram(maximization = False)
    b = p.new_variable(binary = True)
    p.set_objective( sum([b[v] for v in G]) )

    for u in G:
        p.add_constraint( sum([b[v] for v in G[u]]) >= 1 )

    p.solve()
    b = p.get_values(b)
    return len([v for v,i in b.items() if i])


def testiranje_domneve(n):
    grafi = list(graphs.nauty_geng(str(n)+" -c"))
    for graf in grafi:
        if total_domination_number(graf) > (annihilation(graf)) + (powerset_S(graf)):
            show(graf)
            print(annihilation(graf))
            return ("Domneva je ovrzena.")
    return "Domneva ni ovrzena."

#grafi = list(graphs.nauty_geng(str(4)+" -c"))


# GENETIC ALGORITHM

import random
import operator
from math import floor, ceil


# N(t) ~ Poiss(lambda*t)
# N(1) ~ Poiss(lambda)
# N(t) = sum_{k=0...} Ind_{S_k <= t}
# S_k = sum_{i=0...k} T_i; T_i ~ Exp(lambda)
# poissonova porazdelitev iz eksponentne
def poisson(t = 1, lambd = 1/2):
    N = 0
    S = 0
    while S < t:
        N += 1
        S += random.expovariate(lambd)
    return N

# zacetna generacija
# n = stevilo vozlisc
# pop_size = velikost vsake populacije
def generateFirstPopulation(pop_size, n):
    populacija = []
    i = 0
    while i < pop_size:
        osebek = graphs.RandomGNP(n, random.uniform(0, 1))
        if osebek.is_connected():
            populacija.append(osebek)
            i += 1
    return populacija

# hipoteza: tdn(G) <= annihilation(G) + powerset_S(G)
# manjsi fitnes = boljsi graf
def fitness(G):
    return annihilation(G) + powerset_S(G) - total_domination_number(G)

def fitnessPopulation(populacija):
    sez = []
    for osebek in populacija:
        sez.append(fitness(osebek))
    return sez

# uredi populacijo od bolsega grafa do slabsega
# sortirana z elementi (i, j) ... i-ti graf ima fitnes j
def sortPopulation(populacija):
    slovar = {}
    i = 0
    for osebek in populacija:
        slovar[i] = fitness(osebek)
        i += 1
    sortirana = sorted(slovar.items(), key = operator.itemgetter(1))
    nova_populacija = []
    for osebek in sortirana:
        nova_populacija.append(populacija[osebek[0]])
    return nova_populacija

# prvih 'best' je del nove populacije, nekaj 'lucky' je tudi del nove populacije
# best + lucky = velikost populacije
def selectPopulation(populacija, best):
    lucky = pop_size - best
    nova_populacija = []
    sortirana = sortPopulation(populacija)
    for k in range(best):
        nova_populacija.append(sortirana[k])
    srecnezi = random.sample(populacija, k = lucky)
    for osebek in srecnezi:
        nova_populacija.append(osebek)
    random.shuffle(nova_populacija)
    return nova_populacija

# z verjetnostjo 1/3 se doda povezavo, z 1/3 se odstrani in z 1/3 se doda in odstrani
def mutation(osebek):
    prob = random.uniform(0, 1)
    plus_povezave = poisson(lambd = 1/2)
    minus_povezave = poisson(lambd = 1/2)
    if prob <= 1/3:
        for k in range(plus_povezave):
            a, b = osebek.random_vertex(), osebek.random_vertex()
            if a != b:
                osebek.add_edge(a, b)
    elif prob > 1/3 and prob <= 2/3:
        for k in range(minus_povezave):
            a, b = osebek.random_vertex(), osebek.random_vertex()
            osebek.delete_edge(a, b)
            if not osebek.is_connected():
                osebek.add_edge(a, b)
    elif prob > 2/3:
        for k in range(plus_povezave):
            a, b = osebek.random_vertex(), osebek.random_vertex()
            if a != b:
                osebek.add_edge(a, b)
        for k in range(minus_povezave):
            a, b = osebek.random_vertex(), osebek.random_vertex()
            osebek.delete_edge(a, b)
            if not osebek.is_connected():
                osebek.add_edge(a, b)
    return osebek

# z verjetnostjo p izbere posameznika in ga mutira
def mutatePopulation(populacija, p = 5/100):
    nova_populacija = []
    for k in range(len(populacija)):
        rnd = random.random()
        if rnd < p:
            nova_populacija.append(mutation(populacija[k]))
        else:
            nova_populacija.append(populacija[k])
    return nova_populacija

# osebka imata samo enega potomca
# n = stevilo vozlisc
def crossover(n, osebek1, osebek2):
    while True:
        subgraf1 = osebek1.random_subgraph(0.5) #vsebuje neko vozlišče iz osebek1 z verj. 0.5
        subgraf2 = osebek2.random_subgraph(0.5)
        if len(subgraf1.vertices()) + len(subgraf2.vertices()) == n and len(subgraf1.vertices()) >= 1 and len(subgraf1.vertices()) < n and subgraf1.is_connected() and subgraf2.is_connected():
            subgraf1.relabel()
            subgraf2.relabel()
            potomec = subgraf1.disjoint_union(subgraf2)
            # lambda se spreminja z stevilom vozlisc
            nove_povezave = poisson(lambd = log(n/2))
            for k in range(nove_povezave):
                a = subgraf1.random_vertex()
                b = subgraf2.random_vertex()
                potomec.add_edge((0, a), (1, b))
            potomec.relabel()
            break
    return potomec

# stara generacija + potomci
# successful have more chance to mate
def offspringPopulation(populacija):
    #utezi = []
    #M = max(fitnessPopulation(populacija))
    #for k in fitnessPopulation(populacija):
        #utezi.append(M + 1 - k)
    nova_populacija = populacija
    stevilo_parjenj = poisson(t = pop_size)
    for k in range(stevilo_parjenj):
        # starsa = random.choices(populacija, weights = utezi, k = 2) not working ?
        starsa = random.sample(populacija, k = 2)
        nova_populacija.append(crossover(n, starsa[0], starsa[1]))
    return nova_populacija

# n = stevilo vozlisc
# pop_size = velikost populacije
# best = koliko najboljsih izberemo za reprodukcijo
def testGA(n, pop_size, best, stevilo_generacij):
    populacija = generateFirstPopulation(pop_size, n)
    populacija = selectPopulation(populacija, best)
    fitnes = fitnessPopulation(populacija)
    fitnes.sort()
    print(fitnes)
    i = 1
    while i <= stevilo_generacij:
        populacija = offspringPopulation(populacija)
        populacija = mutatePopulation(populacija)
        populacija = selectPopulation(populacija, best)
        fitnes = fitnessPopulation(populacija)
        fitnes.sort()
        print(fitnes)
        if min(fitnes) < 0:
            print("Domneva je ovrzena")
            show(sortPopulation(populacija)[0])
            resitev = sortPopulation(populacija)[0]
            resitev = resitev.to_dictionary()
            return resitev
        i += 1
    return "Domneva ni ovrzena."
︡16ad8a5d-2d0a-4078-b8f5-153a7f2cdfd7︡{"done":true}︡











