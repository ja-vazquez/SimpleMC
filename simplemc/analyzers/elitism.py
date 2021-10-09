try:
    from deap import tools
    from deap import algorithms
except:
    import warnings
    warnings.warn("Please install DEAP library if you want to use ga_deap genetic algorithms.")
    try:
        import sys
        sys.exit("Exit.")
    except:
        pass

import re


def eaSimpleWithElitism(population, toolbox, cxpb, mutpb, ngen, stats=None,
             halloffame=None, outputname='deap_output', verbose=__debug__):
    """This algorithm is similar to DEAP eaSimple() algorithm, with the modification that
    halloffame is used to implement an elitism mechanism. The individuals contained in the
    halloffame are directly injected into the next generation and are not subject to the
    genetic operators of selection, crossover and mutation.
    """
    logbook = tools.Logbook()
    logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])

    # Evaluate the individuals with an invalid fitness
    invalid_ind = [ind for ind in population if not ind.fitness.valid]
    fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
    for ind, fit in zip(invalid_ind, fitnesses):
        ind.fitness.values = fit

    if halloffame is None:
        raise ValueError("halloffame parameter must not be empty!")

    halloffame.update(population)
    hof_size = len(halloffame.items) if halloffame.items else 0

    record = stats.compile(population) if stats else {}
    logbook.record(gen=0, nevals=len(invalid_ind), **record)
    if verbose:
        print(logbook.stream)

    # Write output file on the fly
    f = open('{}_1.txt'.format(outputname), 'w')
    f.write("#Generation(first column) fitness(second column) individual\n")

    # Begin the generational process
    for gen in range(1, ngen + 1):

        # Select the next generation individuals
        offspring = toolbox.select(population, len(population) - hof_size)

        # Vary the pool of individuals
        offspring = algorithms.varAnd(offspring, toolbox, cxpb, mutpb)

        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        # SimpleMC_change
        # Write generation, individual and fitness in output file
        fitnesses_all = toolbox.map(toolbox.evaluate, offspring)
        for indall, fitall in zip(offspring, fitnesses_all):
            strindall = str(indall).lstrip('[').rstrip(']')
            strfitall = str(fitall).lstrip('(').rstrip(')')
            strrow = "{} {} {}\n".format(gen, strfitall, strindall)
            strrow = re.sub(',', '', strrow)
            f.write(strrow)

        # add the best back to population:
        offspring.extend(halloffame.items)

        # Update the hall of fame with the generated individuals
        halloffame.update(offspring)

        # Replace the current population by the offspring
        population[:] = offspring

        # Append the current generation statistics to the logbook
        record = stats.compile(population) if stats else {}
        logbook.record(gen=gen, nevals=len(invalid_ind), **record)
        if verbose:
            print(logbook.stream)

    f.close()
    return population, logbook, gen

