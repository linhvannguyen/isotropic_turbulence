function ids = border_per (ids, minid, maxid) 
%border_per handle boundary by periodicity
ids(ids<minid) = maxid + ids(ids<minid);
ids(ids>maxid) = ids(ids>maxid) - maxid;
end