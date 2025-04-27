#include "algebra.h"
#include <stdio.h>
#include <math.h>

Matrix create_matrix(int row, int col)
{
    Matrix m;
    m.rows = row;
    m.cols = col;
    return m;
}

Matrix add_matrix(Matrix a, Matrix b)
{
    // 检查行列是否匹配
    if (a.rows != b.rows || a.cols != b.cols) {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);  // 返回空矩阵
    }

    // 创建结果矩阵（行数和列数与输入一致）
    Matrix result = create_matrix(a.rows, a.cols);

    // 逐元素相加
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[i][j] = a.data[i][j] + b.data[i][j];
        }
    }

    return result;
}

Matrix sub_matrix(Matrix a, Matrix b)
{
    // 检查行列是否匹配
    if (a.rows != b.rows || a.cols != b.cols) {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);  // 返回空矩阵
    }

    // 创建结果矩阵（行数和列数与输入一致）
    Matrix result = create_matrix(a.rows, a.cols);

    // 逐元素相减
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[i][j] = a.data[i][j] - b.data[i][j];
        }
    }

    return result;
}

Matrix mul_matrix(Matrix a, Matrix b)
{
    // 检查行列是否匹配
    if (a.cols != b.rows) {
        printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n");
        return create_matrix(0, 0);  // 返回空矩阵
    }

    // 创建结果矩阵（行数和列数与输入一致）
    Matrix result = create_matrix(a.rows, b.cols);

    // 逐元素相乘
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < b.cols; j++) {
            for (int k = 0; k < a.cols; k++) {
                result.data[i][j] += a.data[i][k] * b.data[k][j];
            }
        }
    }

    return result;
}

Matrix scale_matrix(Matrix a, double k)
{
    for(int i = 0; i < a.rows; i++) {
        for(int j = 0; j < a.cols; j++) {
            a.data[i][j] *= k;
        }
    }
    return a;
}

Matrix transpose_matrix(Matrix a)
{
    // 创建一个新的矩阵，其行数和列数与原矩阵相反
    Matrix transposed = create_matrix(a.cols, a.rows);

    for(int i = 0; i < a.rows; i++) {
        for(int j = 0; j < a.cols; j++) {
            // 将原矩阵的元素转置到新矩阵中
            transposed.data[j][i] = a.data[i][j];
        }
    }
    return transposed;
}


double det_matrix(Matrix a) {
    // 检查是否为方阵
    if (a.rows != a.cols) {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
    // 基本情况：1x1 矩阵直接返回元素值
    if (a.rows == 1) {
        return a.data[0][0];
    }
    
    double det = 0;
    int sign = 1;  // 符号因子初始化为正
    
    // 按第0行展开
    for (int col = 0; col < a.cols; col++) {
        // 生成余子式矩阵（去掉第0行和第col列）
        Matrix sub = create_matrix(a.rows - 1, a.cols - 1);
        int sub_row = 0;
        
        // 遍历原矩阵的行（从1开始，跳过第0行）
        for (int i = 1; i < a.rows; i++) {
            int sub_col = 0;
            // 遍历原矩阵的列，跳过第col列
            for (int j = 0; j < a.cols; j++) {
                if (j == col) continue;
                sub.data[sub_row][sub_col] = a.data[i][j];
                sub_col++;
            }
            sub_row++;
        }
        
        // 累加：当前元素 * 符号元素 * 余子式的行列式
        det += sign * a.data[0][col] * det_matrix(sub);
        sign *= -1;  // 交替符号
    }
    return det;
}

Matrix inv_matrix(Matrix a) {
    // 检查是否为方阵
    if (a.rows != a.cols) {
        printf("Error: The matrix must be a square matrix.\n");
        return create_matrix(0, 0);
    }
    int n = a.rows;

    // 计算行列式
    double det = det_matrix(a);
    if (det == 0) {
        printf("Error: The matrix is singular.\n");
        return create_matrix(0, 0);
    }

    // 生成伴随矩阵
    Matrix adj = create_matrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // 生成去掉第 j 行和第 i 列的子矩阵
            Matrix sub = create_matrix(n-1, n-1);
            int sub_row = 0;
            for (int row = 0; row < n; row++) {
                if (row == j) continue;  // 跳过第 j 行
                int sub_col = 0;
                for (int col = 0; col < n; col++) {
                    if (col == i) continue;  // 跳过第 i 列
                    sub.data[sub_row][sub_col] = a.data[row][col];
                    sub_col++;
                }
                sub_row++;
            }
            // 计算余子式并乘以符号元素 (-1)^(i+j)
            double cofactor = det_matrix(sub);
            int sign = ((i + j) % 2 == 0) ? 1 : -1;
            adj.data[i][j] = sign * cofactor;
        }
    }

    // 计算逆矩阵：伴随矩阵每个元素除以行列式
    Matrix inv = create_matrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inv.data[i][j] = adj.data[i][j] / det;
        }
    }

    return inv;
}

int rank_matrix(Matrix a) {
    int rank = 0;
    const double EPSILON = 1e-10;  // 浮点精度阈值

    // 遍历所有列，直到处理完所有列或行
    for (int col = 0; col < a.cols && rank < a.rows; col++) {
        // 寻找当前列中绝对值最大的行（从 rank 行开始）
        int max_row = rank;
        double max_val = fabs(a.data[rank][col]);
        for (int i = rank + 1; i < a.rows; i++) {
            if (fabs(a.data[i][col]) > max_val) {
                max_row = i;
                max_val = fabs(a.data[i][col]);
            }
        }

        // 若当前列全为 0，跳过
        if (max_val < EPSILON) {
            continue;
        }

        // 交换当前行与最大行
        if (max_row != rank) {
            for (int k = col; k < a.cols; k++) {
                double temp = a.data[rank][k];
                a.data[rank][k] = a.data[max_row][k];
                a.data[max_row][k] = temp;
            }
        }

        // 消去当前列下方的所有行
        for (int i = rank + 1; i < a.rows; i++) {
            double factor = a.data[i][col] / a.data[rank][col];
            for (int j = col; j < a.cols; j++) {
                a.data[i][j] -= factor * a.data[rank][j];
                // 处理浮点误差，将极小值视为 0
                if (fabs(a.data[i][j]) < EPSILON) {
                    a.data[i][j] = 0.0;
                }
            }
        }

        // 秩增加
        rank++;
    }

    return rank;
}

double trace_matrix(Matrix a) {
    // 检查是否为方阵
    if (a.rows != a.cols) {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
    
    // 计算对角线元素之和
    double trace = 0.0;
    for (int i = 0; i < a.rows; i++) {
        trace += a.data[i][i];
    }
    return trace;
}

void print_matrix(Matrix a)
{
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            // 按行打印，每个元素占8个字符的宽度，小数点后保留2位，左对齐
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}