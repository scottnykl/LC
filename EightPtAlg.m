%X is a set of 8 3x1 vectors of normalized pixels that correspond to a set of 8 3x1 normalized vectors
%in Y. X is from the head image and Y is from the trailing image.
%Assumes: B'*E*A = 0.
function E_p = EightPtAlg( A, B )

%We need to build the Y matrix. Y will be an 8x9 matrix where each column 
%represents the Kth matching points between the head and trail images.

   Y = zeros(9,10);

   for m = 1 : size(Y,2)
      Y(1,m) = B(1,m)*A(1,m);
      Y(2,m) = B(1,m)*A(2,m);
      Y(3,m) = B(1,m);
      Y(4,m) = B(2,m)*A(1,m);
      Y(5,m) = B(2,m)*A(2,m);
      Y(6,m) = B(2,m);
      Y(7,m) = A(1,m);
      Y(8,m) = A(2,m);
      Y(9,m) = -1;
   end

   
   %Now we have Y. So we compute the smallest singular value of Y.
   %Find the "least eigenvector, f of Y'Y.
   A = Y';
      
   %[U,S,V] = svds(Y,1,.000000001); %Computes singular vector closest to .000000001LABEL XYZ
   [U,S,V] = svd(A); %Computes singular vector closest to .000000001
            %Columns of U are eigenvectors of YY'
            %Columns of V are eigenvectors of Y'Y
   
   %E_est = vec2mat( U, 3 ); %may need to TRANSPOSE THIS LABEL XYZ
   E_est = vec2mat( V(:,9), 3 ); %may need to TRANSPOSE THIS, take 9th column
   
   %U' * Y %This shows that the left singular values are prependicular to
   %the matched points.
   
   %Now we need to enforce the internal contstraints which says an
   %Essential matrix has two degrees of freedom -- two singular values are
   %equal and the 3rd is zero.
%We have compute E_est, now we need to do step 3 which involves the Frobenius matrix norm

 [uu,ss,vv] = svd(E_est);
 
 %This diagMat is recommened in the wikipedia 8pt algorithm.
 diagMat = [ ( ss(1,1) + ss(2,2) ) / 2,                         0, 0;
                                     0, ( ss(1,1) + ss(2,2) ) / 2, 0;
                                     0,                         0, 0];
                                 
 %diagMat = [ 1, 0, 0;
 %            0, 1, 0;
 %            0, 0, 0];

 E_p = uu * diagMat * vv';
   
end