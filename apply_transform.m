function [out_x, out_y] = apply_transform(x, y, T)

out_x = (T(1, 1)*x + T(1, 2)*y + T(1, 3)) ./ ...
        (T(3, 1)*x + T(3, 2)*y + T(3, 3));

out_y = (T(2, 1)*x + T(2, 2)*y + T(2, 3)) ./ ...
        (T(3, 1)*x + T(3, 2)*y + T(3, 3));