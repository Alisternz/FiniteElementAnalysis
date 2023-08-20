import numpy as np

def findNextNode(x, y, length, angle):
    angleRad = np.radians(angle)
    nextNodeX = length * np.cos(angleRad) + x
    nextNodeY = length * np.sin(angleRad) + y
    return(nextNodeX, nextNodeY)


def main():
    print(findNextNode(0,0,5,53.13))


main()
