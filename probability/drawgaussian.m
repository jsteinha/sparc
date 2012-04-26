function drawgaussian(mean,cov,color)
    P = chol(cov);
    theta = (0:.001:1)*2*pi;
    pts = P*[cos(theta);sin(theta)];
    plot(pts(1,:)+mean(1),pts(2,:)+mean(2),color);
end