import infection_model
from matplotlib import pyplot as plt

def generate_random_indices(infec_m, num):

    rand_x = []
    rand_y = []

    for i in range(0, round(num / 2) ):
        ind_tup = infection_model.rand_ind1(infec_m)
        rand_x.append(ind_tup[0])
        rand_y.append(ind_tup[1])

    return rand_x, rand_y

def main(times):
    model = infection_model.infection_model(m = 10, n = 10, k = 4, tau = 0.2, delta = 0.1)
    rx, ry = generate_random_indices(model, times)

    plt_name = str(times) + ' Random Indices'

    fig, axs = plt.subplots(1, 2, sharey = True, tight_layout = True)

    axs[0].hist(rx)
    #fig.suptitle(plt_name)
    axs[0].set( xlabel = "Index Value", ylabel = "Frequency")
    axs[0].set_title("Random X Coordinates")

    axs[1].hist(ry)
    axs[1].set( xlabel = "Index Value")
    axs[1].set_title("Random Y Coordinates")
    #axs[0].xlabel("Index Value")

    plt.show()
    

runs = int(input("How many runs: "))
main(runs)
