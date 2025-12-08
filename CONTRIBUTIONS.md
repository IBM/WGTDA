🧱 Extending with New Methods

Add your implementation inside:

src/correlation/<new_method>.py


Then register it in the factory:

elif method == "new_method":
    return compute_new_method(df)