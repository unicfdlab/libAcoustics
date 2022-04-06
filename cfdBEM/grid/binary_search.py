"""Binary search in mesh nodes"""

def binarySearch (array, number, tolerance, first, last):
    """
    Basic binary search recursive algorithm.

    Returns the position of the first encountered number in array which is sufficiently close 
    to target number (according to tolerance) or -1.

    Parameters
    ----------
    array : list of float
    Sorted array
    number : float
        Target number
    tolerance: float
        Tolerance
    first : integer
        Left marker for search
    last : integer
        Right marker for search

    """

    if (first > last):
        return -1
    else:
        midPos = (first + last)//2

        if (abs(number - array[midPos]) > tolerance):
            if (number - array[midPos] < 0):
                return binarySearch (array, number, tolerance, first, midPos - 1)
            else:
                return binarySearch (array, number, tolerance, midPos + 1, last)
        else:
            return midPos

def binarySearchFirst (array, number, pos, tolerance, first, last):
    """
    Find the minimal value in array which is sufficiently close
    to target number.

    Returns the position of encountered value.

    Parameters
    ----------
    array : list of float
    Sorted array
    number : float
        Target number
    pos : integer
    Current position in array
    tolerance: float
        Tolerance
    first : integer
        Left marker for search
    last : integer
        Right marker for search

    """

    if (first >= last or pos == 0 or abs(number - array[pos-1]) > tolerance):
        return pos
    elif (len(array) == 2):
        if abs(number - array[pos-1]) > tolerance:
            return pos
        else:
            return pos - 1
    else:
        midPos = (first + last)//2

        if (abs(number - array[midPos]) > tolerance):
            nextPos = binarySearch (array, number, tolerance, midPos + 1, pos - 1)
            return binarySearchFirst (array, number, nextPos, tolerance, midPos + 1, nextPos)
        else:
            return binarySearchFirst (array, number, midPos, tolerance, first, midPos)

def binarySearchLast (array, number, pos, tolerance, first, last):
    """
    Find the maximal value in array which is sufficiently close
    to target number.

    Returns the position of encountered value.

    Parameters
    ----------
    array : list of float
    Sorted array
    number : float
        Target number
    pos : integer
        Current position in array
    tolerance: float
        Tolerance
    first : integer
        Left marker for search
    last : integer
        Right marker for search

    """

    if (first >= last or pos == len(array)-1 or abs(number - array[pos+1]) > tolerance):
        return pos
    else:
        midPos = (first + last)//2

    if (abs(number - array[midPos]) > tolerance):
        nextPos = binarySearch (array, number, tolerance, pos + 1, midPos - 1)
            return binarySearchLast (array, number, nextPos, tolerance, nextPos, midPos - 1)
    else:
        return binarySearchLast (array, number, midPos+1, tolerance, midPos+1, last)

def coordBinarySearch (sortArrX, coordX, point, tolerance):
    """
    Binary search of duplicate nodes in mesh

    Parameters
    ----------

    sortArrX : tuple of points (number + array of float coordinates)
    Mesh nodes sorted by x-coord
    coordX : list of float
        x-coordinates of points (using just for increasing of computational speed)
    point : array of three float coordinates
    Coordinates of target point
    tolerance :
    Merge tolerance
    """

    if (len(sortArrX) == 0):
        return -1

    else:

        # find x-coord
        xPos = binarySearch(coordX,point[0],tolerance, 0, len(coordX)-1)

        if (xPos == -1):
            return -1
        else:
            xLeft = binarySearchFirst(coordX,point[0],xPos,tolerance, 0, xPos)
            xRight = binarySearchLast(coordX,point[0],xPos,tolerance, xPos, len(coordX)-1)

            # sort by y-coord
            sortArrY = tuple( sorted( sortArrX[xLeft:xRight+1], key = lambda t: t[1][1] ))

            transposeY = zip(*sortArrY)
            coordY = zip(*(transposeY[1]))[1]

            yPos = binarySearch(coordY,point[1],tolerance,0,len(coordY) - 1)

            if (yPos == -1):
                return -1
            else:
                yLeft = binarySearchFirst(coordY,point[1],yPos,tolerance, 0, yPos)
                yRight = binarySearchLast(coordY,point[1],yPos,tolerance, yPos, len(coordY)-1)

                # sort by z-coord
                sortArrZ = tuple( sorted( sortArrY[yLeft:yRight+1], key = lambda t: t[1][2] ))

                transposeZ = zip(*sortArrZ)
                coordZ = zip(*(transposeZ[1]))[2]

                zPos = binarySearch(coordZ,point[2],tolerance,0,len(coordZ)-1)

                if (zPos == -1):
                    return -1
                else:
                    return transposeZ[0][zPos]
