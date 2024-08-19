function X = test_fft(x)
	z = x(1:2:end) + 1j*x(2:2:end);
	Z = fft(z);

	E = Z;
	E(1) = real(E(1));
	E(2:end) += conj(flip(E(2:end)));
	E(2:end) /= 2;

	O = Z;
	O(1) = imag(O(1));
	O(2:end) -= conj(flip(O(2:end)));
	O(2:end) *= -1j/2;

	X = [E + O .* exp(-1j*2*pi*(0:(length(E)-1))/length(x)), E(1) - O(1)];
end
