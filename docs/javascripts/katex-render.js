function renderMath(element) {
  if (typeof renderMathInElement !== "function") {
    return;
  }
  renderMathInElement(element, {
    delimiters: [
      {left: "\\(", right: "\\)", display: false},
      {left: "\\[", right: "\\]", display: true},
      {left: "$$", right: "$$", display: true},
      {left: "$", right: "$", display: false}
    ],
    throwOnError: false,
    errorColor: "#cc0000"
  });
}

document.addEventListener("DOMContentLoaded", () => {
  renderMath(document.body);
});

if (typeof document$ !== "undefined") {
  document$.subscribe(() => {
    renderMath(document.body);
  });
}
