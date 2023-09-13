"""
Module that compares our method to Shrake with a plot
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def comp_plot(dtframe):
    """
    Plot the accessibility values obtained by our program and the ones
    given by Shrake
    Arguments:
        dtframe: A DataFrame containing "shrake_accessible_surface" and "Accessible_surface" columns.
    """
    dtframe["Difference"] = abs(dtframe["shrake_accessible_surface"] - dtframe["Accessible_surface"])
    percentage = 100 * (1 - (dtframe["Difference"].sum() / dtframe["shrake_accessible_surface"].sum()))
    print(f"The percentage of identity between the two methods is : {percentage:.2f}%")
    
    shrake_data = np.array(dtframe["shrake_accessible_surface"])
    our_data = np.array(dtframe["Accessible_surface"])

    plt.plot(shrake_data, label="Shrake result", color="red")
    plt.plot(our_data, label="Our result", color="blue")
    plt.xticks(range(len(shrake_data)), dtframe.index.get_level_values(1))
    plt.xlabel("Residues")
    plt.ylabel("Accessibility")
    plt.legend()
    return plt.savefig("results/plot.png")


if __name__ == "__main__":
    import plot_comparison
    print(help(plot_comparison))

