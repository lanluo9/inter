function res = myReconstruction2(h, X, step, tmpImageSize)
    % preparation. X consists of original image surounded by zeros
    filterSize = size(h);

    numFilter = ceil(((tmpImageSize(1) - filterSize(1))/2 + 1)/step) * 2 + 1;
    imageSize = step * (numFilter - 1) + filterSize(1);
    tmpRes = zeros(imageSize);


    startingPoints(1,:) = 1: step: imageSize - filterSize(1) + 1;
    startingPoints(2,:) = 1: step: imageSize - filterSize(2) + 1;

    for ii = 1: size(startingPoints,2)
        for jj = 1: size(startingPoints,2)

            offset(1) = startingPoints(1,ii);
            offset(2) = startingPoints(2,jj);

            filterLength(1,:) = 0: filterSize(1)-1;
            filterLength(2,:) = 0: filterSize(2)-1;

            Xind = offset(1) + filterLength(1,:);
            Yind = offset(2) + filterLength(2,:);

            tmpRes(Xind, Yind) = tmpRes(Xind, Yind) + X(ii,jj) * h;
        end
    end

    tx = imageSize / 2 - tmpImageSize/2 + 1: imageSize / 2 + tmpImageSize/2;
    res = tmpRes(tx, tx);
