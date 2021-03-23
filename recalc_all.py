from multiprocessing import Process
import os


def casa():
    print 'running features\n\n'
    os.system('python features.py data')
    print 'running casa\n\n'
    os.system('python casa.py no_plots')
    
def individual():
    print 'running individual\n\n'
    os.system('python individual.py data')
    print 'running individual_velocity'
    os.system('python individual_velocity.py no_plots')

if __name__ == '__main__':
    print 'running window_functions\n\n'
    os.system('python window_functions.py data')
    
    casa_p=Process(target=casa)
    ind_p=Process(target=individual)
    
    casa_p.start()
    ind_p.start()
    ind_p.join()
    casa_p.join()
    
    print 'running table_maker\n\n'
    os.system('python table_maker.py')
    print 'running deluxe_table_maker\n\n'
    os.system('python deluxe_table_maker.py')
    print 'running make_plots\n\n'
    os.system('python make_plots.py')

    
print 'done'