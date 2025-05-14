function angle = vecAngle(u,v)
    cosAngle = dot(u,v)/(norm(u)*norm(v));
    angle = acos(cosAngle);
end