function foo(temp)
	minval=temp;
	maxval = temp;
	for i=1:ndims(temp)
		minval = min(minval);
		maxval = max(maxval);
	end
	disp(['min: ' num2str(minval) ', max: ' num2str(maxval)]);
end
