//! Declarative `#[pyclass]` inner wrapper boilerplate.

/// `{ inner: T }` pyclass payload + `wrap(inner)`.
#[macro_export]
macro_rules! pyhandle {
    (
        $(#[$meta:meta])*
        $vis:vis struct $pyname:ident {
            inner: $inner:ty $(, $($rest:tt)*)?
        }
    ) => {
        $(#[$meta])*
        $vis struct $pyname {
            inner: $inner,
        }

        impl $pyname {
            pub fn wrap(inner: $inner) -> Self {
                Self { inner }
            }

            pub fn borrow_inner(&self) -> &$inner {
                &self.inner
            }
        }
    };
}
