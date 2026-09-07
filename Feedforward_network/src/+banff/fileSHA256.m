function value = fileSHA256(path)
%FILESHA256 Bind a saved dataset model to the actual prepared data file.
fid = fopen(path,'rb');
assert(fid >= 0,'Cannot read dataset: %s',path);
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
digest = java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid)
    chunk = fread(fid,1024*1024,'*uint8');
    digest.update(typecast(chunk,'int8'));
end
bytes = typecast(digest.digest(),'uint8');
value = lower(reshape(dec2hex(bytes,2).',1,[]));
end
