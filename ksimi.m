function kernelsimilarity = ksimi(f,kernel)  
    fnorm=f/norm(f,'fro');
    knorm=kernel/norm(kernel,'fro');
    kernelsimilarity=fnorm(:)'*knorm(:);