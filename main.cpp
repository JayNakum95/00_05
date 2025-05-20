#include "Novice.h"
#include <cmath>
#include <cstdio>

const char kWindowTitle[] = "GC2B_07_ナクム_ジェイ_ハルシュバルダン";

// 行列・ベクトル定義
struct Matrix4x4 {
    float m[4][4];
};
struct Vector3 {
    float x, y, z;
};
// 画面表示用定数
static const int kRowHeight = 20;
static const int kColumnWidth = 60;

// 行列の加算
static Matrix4x4 Add(Matrix4x4& m1, Matrix4x4& m2) {
    Matrix4x4 result{};
    for (int row = 0; row < 4; ++row)
        for (int col = 0; col < 4; ++col)
            result.m[row][col] = m1.m[row][col] + m2.m[row][col];
    return result;
}

// 行列の減算
static Matrix4x4 Subtract(Matrix4x4& m1, Matrix4x4& m2) {
    Matrix4x4 result{};
    for (int row = 0; row < 4; ++row)
        for (int col = 0; col < 4; ++col)
            result.m[row][col] = m1.m[row][col] - m2.m[row][col];
    return result;
}

// 行列の乗算
static Matrix4x4 Multiply(Matrix4x4& m1, Matrix4x4& m2) {
    Matrix4x4 result{};
    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 4; ++col) {
            result.m[row][col] = 0;
            for (int k = 0; k < 4; ++k) {
                result.m[row][col] += m1.m[row][k] * m2.m[k][col];
            }
        }
    }
    return result;
}


// 行列の転置
static Matrix4x4 Transpose(const Matrix4x4& m) {
    Matrix4x4 result{};
    for (int row = 0; row < 4; ++row)
        for (int col = 0; col < 4; ++col)
            result.m[row][col] = m.m[col][row];
    return result;
}

// 単位行列
static Matrix4x4 MakeIdentity4x4() {
    Matrix4x4 result{};
    for (int i = 0; i < 4; ++i)
        result.m[i][i] = 1.0f;
    return result;
}

static Matrix4x4 Inverse(Matrix4x4& m) {
    Matrix4x4 result = MakeIdentity4x4();
    Matrix4x4 temp = m;

    for (int i = 0; i < 4; ++i) {
        if (fabsf(temp.m[i][i]) < 1e-6f) {
            bool swapped = false;
            for (int j = i + 1; j < 4; ++j) {
                if (fabsf(temp.m[j][i]) > 1e-6f) {
                    for (int k = 0; k < 4; ++k) {
                        std::swap(temp.m[i][k], temp.m[j][k]);
                        std::swap(result.m[i][k], result.m[j][k]);
                    }
                    swapped = true;
                    break;
                }
            }
            if (!swapped) {
                // 逆行列が存在しない
                return MakeIdentity4x4(); // 代替措置
            }
        }

        // 対角要素を1にする
        float diag = temp.m[i][i];
        for (int k = 0; k < 4; ++k) {
            temp.m[i][k] /= diag;
            result.m[i][k] /= diag;
        }

        // 他の行のi列を0にする
        for (int j = 0; j < 4; ++j) {
            if (i == j) continue;
            float factor = temp.m[j][i];
            for (int k = 0; k < 4; ++k) {
                temp.m[j][k] -= factor * temp.m[i][k];
                result.m[j][k] -= factor * result.m[i][k];
            }
        }
    }

    return result;
}

Matrix4x4 MakeRotationX(float radian) {
    Matrix4x4 result = MakeIdentity4x4();
    float cosradian = cosf(radian);
    float sinradian = sinf(radian);
    result.m[1][1] = cosradian;
    result.m[1][2] = sinradian;
    result.m[2][1] = -sinradian;
    result.m[2][2] = cosradian;
    return result;
}
Matrix4x4 MakeRotationY(float radian) {
    Matrix4x4 result = MakeIdentity4x4();
    float cosradian = cosf(radian);
    float sinradian = sinf(radian);
    result.m[0][0] = cosradian;
    result.m[0][2] = -sinradian;
    result.m[2][0] = sinradian;
    result.m[2][2] = cosradian;
    return result;
}
Matrix4x4 MakeRotationZ(float radian) {
    Matrix4x4 result = MakeIdentity4x4();
    float cosradian = cosf(radian);
    float sinradian = sinf(radian);
    result.m[0][0] = cosradian;
    result.m[0][1] = sinradian;
    result.m[1][0] = -sinradian;
    result.m[1][1] = cosradian;
    return result;
}
static Matrix4x4 MakeAffineMatrix(const Vector3& scale, const Vector3& rotate, const Vector3& translate) {
    Matrix4x4 result = MakeIdentity4x4(); 
    Matrix4x4 rotationX = MakeRotationX(rotate.x);
    Matrix4x4 rotationY = MakeRotationY(rotate.y);
    Matrix4x4 rotationZ = MakeRotationZ(rotate.z);


	Matrix4x4 temporaryMatrix = Multiply(rotationY, rotationZ);
    Matrix4x4 rotation = Multiply(rotationX, temporaryMatrix);

    // Apply scaling after rotation
    result.m[0][0] = rotation.m[0][0] * scale.x;
    result.m[0][1] = rotation.m[0][1] * scale.x;
    result.m[0][2] = rotation.m[0][2] * scale.x;

    result.m[1][0] = rotation.m[1][0] * scale.y;
    result.m[1][1] = rotation.m[1][1] * scale.y;
    result.m[1][2] = rotation.m[1][2] * scale.y;

    result.m[2][0] = rotation.m[2][0] * scale.z;
    result.m[2][1] = rotation.m[2][1] * scale.z;
    result.m[2][2] = rotation.m[2][2] * scale.z;

    // Apply translation last
    result.m[3][0] = translate.x;
    result.m[3][1] = translate.y;
    result.m[3][2] = translate.z;

    return result;
}


// 画面に行列を表示
static void MatrixScreenPrintf(int x, int y, Matrix4x4& matrix, const char* label) {
    Novice::ScreenPrintf(x, y - 20, "%s", label);
    for (int row = 0; row < 4; ++row) {
        for (int column = 0; column < 4; ++column) {
            Novice::ScreenPrintf(
                x + column * kColumnWidth,
                y + row * kRowHeight,
                "%6.02f",
                matrix.m[row][column]);
        }
    }
}

// メイン
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {
    Novice::Initialize(kWindowTitle, 1280, 720);

    Matrix4x4 m1 = {
        3.2f, 0.7f, 9.6f, 4.4f,
        5.5f, 1.3f, 7.8f, 2.1f,
        6.9f, 8.0f, 2.6f, 1.0f,
        0.5f, 7.2f, 5.1f, 3.3f
    };

    Matrix4x4 m2 = {
        4.1f, 6.5f, 3.3f, 2.2f,
        8.8f, 0.6f, 9.9f, 7.7f,
        1.1f, 5.5f, 6.6f, 0.0f,
        3.3f, 9.9f, 8.8f, 2.2f
    };

    // 計算
    Matrix4x4 resultAdd = Add(m1, m2);
    Matrix4x4 resultSubtract = Subtract(m1, m2);
    Matrix4x4 resultMultiply = Multiply(m1, m2);
    Matrix4x4 inverseM1 = Inverse(m1);  // 仮の単位行列
    Matrix4x4 inverseM2 = Inverse(m2);  // 仮の単位行列
    Matrix4x4 transposeM1 = Transpose(m1);
    Matrix4x4 transposeM2 = Transpose(m2);
    Matrix4x4 identity = MakeIdentity4x4();
	Vector3 scale = { 1.2f, 0.79f, -2.1f };
    Vector3 rotate = { 0.4f,1.43f,-0.8f };
	Vector3 translate = { 2.7f, -4.15f, 1.53f };
    

    //Matrix4x4 rotation = Multiply(rotationX, Multiply(rotationY, rotationZ));
	Matrix4x4 worldMatrix = MakeAffineMatrix(scale, rotate, translate);
    // ループ
    while (Novice::ProcessMessage() == 0) {
        Novice::BeginFrame();

		MatrixScreenPrintf(0, 30, worldMatrix, "worldMatrix");


        Novice::EndFrame();
    }

    Novice::Finalize();
    return 0;
}
