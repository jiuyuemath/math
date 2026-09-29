#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

const double EPS = 1e-10;

void printMatrix(const vector<vector<double>> &matrix)
{
    for (const auto &row : matrix)
    {
        for (size_t j = 0; j < row.size(); ++j)
        {
            cout << (fabs(row[j]) < EPS ? 0 : row[j]);
            if (j + 1 < row.size())
                cout << ' ';
        }
        cout << '\n';
    }
}

int main()
{
    int rows, cols;
    cout << "请输入矩阵的行数和列数：";
    if (!(cin >> rows >> cols) || rows <= 0 || cols <= 0)
    {
        cerr << "行数和列数必须为正整数。\n";
        return 1;
    }

    vector<vector<double>> a(rows, vector<double>(cols));
    cout << "请输入矩阵元素（按行输入）：\n";
    for (auto &row : a)
    {
        for (double &value : row)
        {
            if (!(cin >> value))
            {
                cerr << "矩阵元素输入无效。\n";
                return 1;
            }
        }
    }

    // Gauss-Jordan elimination: convert the input matrix to RREF.
    int pivotRow = 0;
    for (int col = 0; col < cols && pivotRow < rows; ++col)
    {
        int best = pivotRow;
        for (int r = pivotRow + 1; r < rows; ++r)
            if (fabs(a[r][col]) > fabs(a[best][col]))
                best = r;
        if (fabs(a[best][col]) < EPS)
            continue;

        swap(a[pivotRow], a[best]);
        double divisor = a[pivotRow][col];
        for (double &value : a[pivotRow])
            value /= divisor;
        for (int r = 0; r < rows; ++r)
        {
            if (r == pivotRow)
                continue;
            double factor = a[r][col];
            for (int c = 0; c < cols; ++c)
                a[r][c] -= factor * a[pivotRow][c];
        }
        ++pivotRow;
    }
    for (auto &row : a)
        for (double &value : row)
            if (fabs(value) < EPS)
                value = 0;

    cout << "简化阶梯形矩阵：\n";
    printMatrix(a);

    string answer;
    cout << "是否要解方程？（yes/no）：";
    cin >> answer;
    if (answer == "no" || answer == "NO" || answer == "No")
        return 0;
    if (answer != "yes" && answer != "YES" && answer != "Yes")
    {
        cout << "输入无效，程序结束。\n";
        return 0;
    }

    // For solving, treat the last column as b in the augmented matrix [A|b].
    if (cols < 2)
    {
        cout << "增广矩阵至少需要两列（系数列和常数列）。\n";
        return 0;
    }
    const int variables = cols - 1;
    vector<int> pivotColumn(rows, -1);
    for (int r = 0; r < rows; ++r)
    {
        for (int c = 0; c < variables; ++c)
        {
            if (fabs(a[r][c]) > EPS)
            {
                pivotColumn[r] = c;
                break;
            }
        }
        if (pivotColumn[r] == -1 && fabs(a[r][variables]) > EPS)
        {
            cout << "方程组无解。\n";
            return 0;
        }
    }

    vector<bool> isPivot(variables, false);
    for (int c : pivotColumn)
        if (c >= 0)
            isPivot[c] = true;
    bool unique = true;
    for (bool pivot : isPivot)
        if (!pivot)
            unique = false;

    if (unique)
    {
        cout << "方程组的唯一解：\n";
        for (int r = 0; r < rows; ++r)
            if (pivotColumn[r] >= 0)
                cout << "x" << pivotColumn[r] + 1 << " = " << a[r][variables] << '\n';
    }
    else
    {
        cout << "方程组有无穷多解。自由变量为：";
        bool first = true;
        for (int c = 0; c < variables; ++c)
        {
            if (!isPivot[c])
            {
                if (!first)
                    cout << ", ";
                cout << "x" << c + 1;
                first = false;
            }
        }
        cout << '\n';
    }
    return 0;
}