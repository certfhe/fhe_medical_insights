#include "utils-matrices.h"
#include "utils-ptxt.h"
#include <cassert>



Ciphertext<DCRTPoly> maskRow(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    const size_t rowIndex
)
{
    assert(rowIndex >= 0 && rowIndex < matrixSize && "Invalid row index");

    std::vector<double> mask(matrixSize * matrixSize, 0.0);
    for (size_t i = 0; i < matrixSize; i++)
        mask[matrixSize * rowIndex + i] = 1.0;
    
    return c * mask;
}


Ciphertext<DCRTPoly> maskColumn(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    const size_t columnIndex
)
{
    assert(columnIndex >= 0 && columnIndex < matrixSize && "Invalid column index");

    std::vector<double> mask(matrixSize * matrixSize, 0.0);
    for (size_t i = 0; i < matrixSize; i++)
        mask[matrixSize * i + columnIndex] = 1.0;
    
    return c * mask;
}


Ciphertext<DCRTPoly> replicateRow(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize
)
{
    for (size_t i = 0; i < LOG2(matrixSize); i++)
        c += c >> (1 << (LOG2(matrixSize) + i));

    return c;
}


Ciphertext<DCRTPoly> replicateColumn(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize
)
{
    for (size_t i = 0; i < LOG2(matrixSize); i++)
        c += c >> (1 << i);

    return c;
}


Ciphertext<DCRTPoly> sumRows(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    bool maskOutput,
    const size_t outputRow
)
{
    for (size_t i = 0; i < LOG2(matrixSize); i++)
        c += c >> (1 << (LOG2(matrixSize) + i));

    if (maskOutput)
        c = maskRow(c, matrixSize, outputRow);

    return c;
}


Ciphertext<DCRTPoly> sumColumns(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    bool maskOutput
)
{
    for (size_t i = 0; i < LOG2(matrixSize); i++)
        c += c << (1 << i);
    
    if (maskOutput)
        c = maskColumn(c, matrixSize, 0);
    
    return c;
}


Ciphertext<DCRTPoly> transposeRow(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    bool maskOutput
)
{
    for (size_t i = 1; i <= LOG2(matrixSize); i++)
        c += c >> (matrixSize * (matrixSize - 1) / (1 << i));
    
    if (maskOutput)
        c = maskColumn(c, matrixSize, 0);

    return c;
}


Ciphertext<DCRTPoly> transposeColumn(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    bool maskOutput
)
{
    for (size_t i = 1; i <= LOG2(matrixSize); i++)
        c += c << (matrixSize * (matrixSize - 1) / (1 << i));

    if (maskOutput)
        c = maskRow(c, matrixSize, 0);

    return c;
}


std::vector<int32_t> getRotationIndices
(
    const size_t matrixSize
)
{
    std::vector<int32_t> indices;

    int32_t index;
    for (size_t i = 0; i < LOG2(matrixSize); i++)
    {
        index = 1 << i;
        indices.push_back(index);   // sumColumns
        indices.push_back(-index);  // replicateColumn

        index = 1 << (LOG2(matrixSize) + i);
        indices.push_back(-index);  // replicateRow, sumRows

        index = (matrixSize * (matrixSize - 1) / (1 << (i + 1)));
        indices.push_back(index);   // transposeColumn
        indices.push_back(-index);  // transposeRow
    }
     indices.push_back(1);
    return indices;
}


