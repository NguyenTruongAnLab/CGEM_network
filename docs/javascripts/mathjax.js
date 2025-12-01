window.MathJax = {
  tex: {
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex"
  }
};

(function () {
  const typeset = () => {
    MathJax.typesetPromise().catch((error) => {
      console.warn("MathJax typeset failed:", error);
    });
  };

  const attachListeners = () => {
    if (typeof document$ !== "undefined") {
      document$.subscribe(typeset);
    } else {
      window.addEventListener("DOMContentLoaded", typeset);
    }
  };

  MathJax.startup.promise.then(() => {
    typeset();
    attachListeners();
  });
})();
