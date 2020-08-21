def processSTR(s):
    rs = []
    for s in s.split(' '):
        if '.' in s:
            rs.append(float(s))
        else:
            rs.append(int(s))
    return rs

import math
class Matrix(object):
    def __init__(self, array):
        self.array = array
        self.rows, self.cols = 0, 0
        if isinstance(array[0], list):
            self.rows = len(array)
            self.cols = len(array[0])
        else:
            self.rows = 1
            self.cols = len(array)

    @staticmethod
    def delrowcol(arr, r, c): #r, c will be deleted
        row, col = len(arr), len(arr[0])
        res = []
        for i in range(0, row):
            if i != r:
                R = []
                for j in range(0, col):
                    if j!=c:
                        R.append(arr[i][j])
                res.append(R)
        return res

    def __str__0(self):
        s = ""
        for rows in self.array:
            for col in rows:
                s+=str(col)+' '
            s+= '\n'
        return s

    def __str__(self):
        string = ''
        for row in self.array:
            for element in row:
                if element >= 0:
                    addition = str(element).split('.')
                    if len(addition) == 2:
                        addition = f"{addition[0]}.{addition[1][:2]} "
                    else:
                        addition = f"{addition[0]} "
                    string += addition
                else:
                    addition = str(element).split('.')
                    if len(addition) == 2:
                        addition = f"{addition[0]}.{addition[1][:2]} "
                    else:
                        addition = f"{addition[0]} "
                    string += addition
            string += '\n'
        return string

    def __mul__(self, Z):
        if type(Z) == int or type(Z) == float:
            return self.scalar(Z, mult = True)
        #Format A : r * c, B : r * c
        if self.cols != Z.rows:
            return 'fail'
        res = Matrix([[0] * Z.cols for _ in range(self.rows)])
        for L0 in range(res.rows):
            for L1 in range(res.cols):
                for L2 in range(Z.rows):
                    res.array[L0][L1] += self.array[L0][L2] * Z.array[L2][L1]
        return res

    def __add__(self, Z):
        x, y = Z.rows, Z.cols
        res = [[0]*y for _ in range(x)]
        for xi in range(x):
            for yi in range(y):
                res[xi][yi] = self.array[xi][yi] + Z.array[xi][yi]

        return Matrix(res)

    @staticmethod
    def LaplaceDet(arr, rowSelect=0):
        assert len(arr) == len(arr[0]) , 'NonSquareMatrixError'
        if len(arr) == 2 & len(arr[0]) == 2:
            return (arr[0][0] * arr[1][1]) - (arr[1][0] * arr[0][1])
        else:
            det = 0
            sign = -1 if rowSelect % 2 else 1
            for minor in range(len(arr)):
                det += arr[rowSelect][minor] * Matrix.LaplaceDet(Matrix.delrowcol(arr, rowSelect, minor)) * sign
                sign *= -1
            return det

    @classmethod
    def inputMatrix(cls, ordinate):
        ScanRC = input('Enter size of %smatrix: > '%(ordinate)).split(" ")
        print('Enter %smatrix:'%(ordinate))
        Arr = [processSTR(input('> ')) for _ in range(int(ScanRC[0]))]
        return Matrix(Arr)

    @staticmethod
    def applyCofactorSigns(arr):
        for i in range(len(arr)):
            for j in range(len(arr[0])):
                if (i + j) & 1:
                    arr[i][j] *= -1
        
    def transpose(self, TYPE = 1):
        if TYPE == 1:
            x, y = self.rows, self.cols
            res = [[0 for _ in range(x)] for _ in range(y)]
            for i in range(x):
                for j in range(y):
                    res[j][i] = self.array[i][j] 
        
            return Matrix(res)

        if TYPE == 2:
            x, y = self.rows, self.cols
            res = [[0 for _ in range(y)] for _ in range(x)]
            for i in range(x - 1):
                for j in range(y - i - 1):
                    '''
                    temp = self.array[i][j]
                    self.array[i][j] = self.array[x - j - 1][y - i - 1]
                    self.array[x - j - 1][y - i - 1] = temp
                    '''
                    res[x - j - 1][y - i - 1] = self.array[i][j]
                    res[i][j] = self.array[x - j - 1][y - i - 1]
            
            for k in range(x):
                res[k][x - k - 1] = self.array[k][x - k - 1]
                
            return Matrix(res)

        if TYPE == 3:
            x, y = self.rows, self.cols
            res = [[0 for _ in range(y)] for _ in range(x)]
            y_i = y //  2
            for i in range(x):
                for j in range(y_i):
                    res[i][j] = self.array[i][y - j - 1]
                    res[i][y - j - 1] = self.array[i][j]

            return Matrix(res)

        if TYPE == 4:
            x, y = self.rows, self.cols
            res = [[0 for _ in range(y)] for _ in range(x)]
            x_i = x //  2
            for i in range(x_i):
                for j in range(y):
                    res[i][j] = self.array[x - 1 - i][j]
                    res[x - 1 - i][j] = self.array[i][j]

            return Matrix(res)

    def scalar(self, num, mult = False):
        for r in range(self.rows):
            for c in range(self.cols):
                if mult:
                    self.array[r][c] *= num
                else:
                    self.array[r][c] += num

        return self

    def determinant_iter(self):
        row, col = self.rows, self.cols
        assert row == col, 'Non Square Error'
        det = 0
        for core in range(row):
            res, neg_res = 1, 1
            for shift in range(col):
                neg_res *= self.array[row - shift - 1][(core + shift) % row]
                res *= self.array[shift][(core + shift) % row]
            det = det + res - neg_res
        return det

    def inverseMat(self):
        mmR, mmC = self.rows, self.cols
        assert mmR == mmC, 'NonSqareMatrixError'
        mmX = Matrix([[0]*mmC for _ in range(mmR)])
        for i in range(mmR):
            for j in range(mmC):
                mmX.array[i][j] = Matrix.LaplaceDet(Matrix.delrowcol(self.array,i, j))

        Matrix.applyCofactorSigns(mmX.array)
        mmX = mmX.transpose(TYPE=1)
        mmX.scalar(1 / Matrix.LaplaceDet(self.array) , mult=True)
        return mmX

def main():
    while(1):
        p = 1
        print("1. Add matrices\n2. Multiply matrix by a constant\n3. Multiply matrices\n4. Transpose Matrix\n5. Calculate a determinant\n6. Inverse matrix\n0. Exit")
        choice = int(input("Your choice: "))
        f,s = 'first ', 'second '
        res = ''
        if choice == 0:
            break   
        elif(choice == 1): #Add
            A = Matrix.inputMatrix(f)
            B = Matrix.inputMatrix(s)
            res = A+B
        elif(choice == 2): #c
            C = Matrix.inputMatrix(ordinate='')
            C0 = int(input('Enter constant: '))
            res = C * C0
        elif(choice == 3): #c
            C = Matrix.inputMatrix(f)
            D = Matrix.inputMatrix(s)
            res = C * D
        elif(choice==4):
            print('1. Main diagonal\n2. Side diagonal\n3. Vertical line\n4. Horizontal line')
            c = int(input('Your choice: '))
            M = Matrix.inputMatrix(ordinate='')
            res = M.transpose(TYPE=c)
        elif(choice == 5):
            D = Matrix.inputMatrix(ordinate='')
            if D.cols == D.rows and D.cols == 1:
                res = D.array[0][0]
            else:
                res = Matrix.LaplaceDet(D.array)
        elif(choice == 6):
            I = Matrix.inputMatrix(ordinate='')
            det = Matrix.LaplaceDet(I.array)
            if det != 0:
                res = I.inverseMat()
            else:
                res = "This matrix doesn't have an inverse."
                p = 0
        if p:
            print('The result is')
        print(res)
        

if __name__ == "__main__":
    main()