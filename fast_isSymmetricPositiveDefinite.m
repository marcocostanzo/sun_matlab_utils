function isSymPD = fast_isSymmetricPositiveDefinite(A)
    % fast check if a matrix is symmetric positive definite

    try 
        chol(A);
        isSymPD = true;
    catch ME
        isSymPD = false;
    end

end

