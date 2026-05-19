# Step Card

Creates a styled card for an analysis step with a header, body, and an
optional "next step" hint footer.

## Usage

``` r
step_card(title, ..., id = NULL, full_screen = FALSE, next_step = NULL)
```

## Arguments

- title:

  Character. Title for the step.

- ...:

  UI elements for the card body.

- id:

  Optional card ID.

- full_screen:

  Logical. Allow full-screen viewing.

- next_step:

  Optional character. If supplied, a muted footer is appended with a
  right-pointing arrow and this text to guide the user toward the next
  action (e.g. "Proceed to Step 3: Normalization").

## Value

A bslib card.
